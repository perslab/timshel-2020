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


dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/"
setwd(wd)


# ======================================================================= #
# ================================ FUNCTIONS ================================= #
# ======================================================================= #

load_ldsc_cts_results <- function(file.ldsc_cts, dataset_prefix) {
  # dataset_prefix: character vector given the data set prefix. This is used to seperate the sem_name from the annotation_name.
  
  df.ldsc_cts <- read_tsv(file.ldsc_cts) # The last column gives a P-value from a one-sided test that the coefficient is greater than zero. 
  ### Name column examples in CTS file:
  # celltypes.mousebrain.all.binary__mousebrain_all.DEGLU5.tstat
  # celltypes.campbell_lvl1.all__campbell_lvl1.a18.Neurons6.tstat
  # celltypes.tabula_muris.all.binary__tabula_muris.Brain_Non-Myeloid.Bergmann_glial_cell.specificity
  pattern <- sprintf(".*__%s\\.(.*)\\.(.*?)$",dataset_prefix)
  mat.match <- str_match(df.ldsc_cts$Name, pattern) # returns matrix where the first column is the complete match, followed by one column for each capture group
  mat.match
  df.ldsc_cts <- df.ldsc_cts %>% mutate(sem=mat.match[,3],
                                        annotation=mat.match[,2],
                                        dataset=dataset_prefix) %>%
    select(-Name) # drop Name col now
  
  ### split string - ONLY works if there are no dots (.) in the annotation_name. (Works for mousebrain)
  # tmp_split <- stringr::str_split_fixed(df.ldsc_cts$Name,pattern="\\.", n=Inf)
  # df.ldsc_cts <- df.ldsc_cts %>% mutate(sem=tmp_split[,length(tmp_split[1,])],
  #                                       annotation=tmp_split[,length(tmp_split[1,])-1],
  #                                       dataset=tmp_split[,length(tmp_split[1,])-2]) %>% 
  #   select(-Name) # drop Name col now
  
  df.ldsc_cts <- df.ldsc_cts %>% rename(estimate=Coefficient, std.error=Coefficient_std_error, p.value=Coefficient_P_value)
  df.ldsc_cts.tmp.summary <- df.ldsc_cts %>% group_by(sem) %>% summarise(n_obs_sem=n())
  df.ldsc_cts <- left_join(df.ldsc_cts, df.ldsc_cts.tmp.summary, by="sem")
  df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n_obs_sem, true=T, false=F),
                                        p.value.adj = p.value*n_obs_sem)
  # df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
  return(df.ldsc_cts)
}


# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #


### PARAMS
dataset_prefix <- "mousebrain_all"
# dataset_prefix <- "tabula_muris"

# ======================================================================= #
# =============================== LOAD LDSC CTS RESULTS ================================= #
# ======================================================================= #

if (dataset_prefix == "campbell_lvl1") {
  genomic_annotation_perfix <- "celltypes.campbell_lvl1.all"
} else if (dataset_prefix == "campbell_lvl2") {
  genomic_annotation_perfix <- "celltypes.campbell_lvl2.all"
} else if (dataset_prefix == "mousebrain_all") {
  genomic_annotation_perfix <- "celltypes.mousebrain.all"
} else if (dataset_prefix == "tabula_muris") {
  genomic_annotation_perfix <- "celltypes.tabula_muris.all"
} else {
  stop("wrong dataset_prefix")
}


### Load - CONDITIONAL
dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern=sprintf("%s__BMI_Yengo2018(.*).cell_type_results.txt", genomic_annotation_perfix))
# filenames <- filenames[grepl(pattern="__CONDITIONAL__", filenames, perl=T)] # exclude any conditional results
list.dfs <- lapply(file.path(dir.data, filenames), load_ldsc_cts_results, dataset_prefix)
# stringr::str_match(filenames, pattern=sprintf("%s__(.*?)(__CONDITIONAL__)?(.*?).cell_type_results.txt", genomic_annotation_perfix))
tmp <- stringr::str_match(filenames, pattern=sprintf("%s__(.*?)(__CONDITIONAL__)(.*?).cell_type_results.txt", genomic_annotation_perfix))[,4]
tmp[is.na(tmp)] <- "baseline"
tmp
names(list.dfs) <- tmp
names(list.dfs)
df.ldsc_cts <- list.dfs %>% bind_rows(.id="gwas")
df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")

# df.ldsc_cts <- df.ldsc_cts %>% select(gwas, p.value, annotation) %>% spread(key=gwas, value=p.value)

# ======================================================================= #
# ============================== READ METADATA =============================== #
# ======================================================================= #

if (dataset_prefix == "campbell_lvl1") {
  file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/hypothalamus_campbell/campbell_lvl1.cell_type_metadata.csv"
  df.metadata <- read_csv(file.metadata) %>% rename(annotation = cell_type_all_lvl1)
  df.metadata <- df.metadata %>% mutate(color_by_variable = taxonomy_lvl1)
} else if (dataset_prefix == "campbell_lvl2") {
  file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/hypothalamus_campbell/campbell_lvl2.cell_type_metadata.csv"
  df.metadata <- read_csv(file.metadata) %>% rename(annotation = cell_type_all_lvl2)
  df.metadata <- df.metadata %>% mutate(color_by_variable = taxonomy_lvl1)
} else if (dataset_prefix == "mousebrain_all") {
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
  df.metadata <- df.metadata %>% mutate(color_by_variable = Class)
} else if (dataset_prefix == "tabula_muris") {
  file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/tabula_muris/tabula_muris_facs.tissue_celltype.celltype_metadata.csv"
  df.metadata <- read_csv(file.metadata)
  df.metadata <- df.metadata %>% mutate(
    annotation = tissue_celltype,
    color_by_variable = tissue)
  df.metadata <- df.metadata %>% mutate(annotation = stringr::str_replace_all(annotation, pattern="\\s+", replacement="_"))
} else {
  stop("wrong dataset_prefix")
}

df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data


# ======================================================================= #
# =============================== EXPORT ================================= #
# ======================================================================= #

### Multi CONDITION
# df.ldsc_cts %>% write_csv(sprintf("out.tmp.190111.%s.sem_mean.CONDITIONAL_BMI_Yengo2018.csv", dataset_prefix))



# ======================================================================= #
# ================================ PLOT: cell priori [FACET WRAP MULTI-CONDITIONAL] ================================= #
# ======================================================================= #


df.plot <- df.ldsc_cts
fdr_threshold <- 0.05/nrow(df.plot %>% filter(gwas==(df.plot %>% pull(gwas))[1])) # *OBS*: TMP DIRTY CODE FOR CONDITIONAL
p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
  geom_point(aes(color=color_by_variable, size=estimate)) + 
  ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
  geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
  facet_wrap(~gwas, ncol = 1) + 
  labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
  # ggtitle(sprintf("%s - %s", "BMI_Yengo2018", name_elem)) +
  # ggtitle(sprintf("%s - %s", dataset_prefix, name_elem)) +
  theme_classic() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p
file.out <- sprintf("out.tmp.190111.plot.cell_prioritization.%s.%s.sem_mean.color_by_class.pdf", dataset_prefix, "CONDITIONAL")
ggsave(p, filename=file.out, width=15, height=20)



# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #



# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #




# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

