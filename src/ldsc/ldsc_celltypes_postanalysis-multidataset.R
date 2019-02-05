############### SYNOPSIS ###################
# Analyze and export LDSC results for multiple datasets (tabula muris, mousebrain, ...)


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
# ============================ PARAMETERS =============================== #
# ======================================================================= #


### PARAMS
dataset_prefix <- "mousebrain_all"
# dataset_prefix <- "tabula_muris"
# dataset_prefix <- "campbell_lvl1"
# dataset_prefix <- "campbell_lvl2"

genomic_annotation_prefix <- get_genomic_annotation_prefix(dataset_prefix)

# ======================================================================= #
# =============================== LOAD LDSC CTS RESULTS ================================= #
# ======================================================================= #

### Single GWAS
file.ldsc_cts <- sprintf("/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__BMI_UPDATE_Yengo2018.cell_type_results.txt", genomic_annotation_prefix)
# file.ldsc_cts <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.mousebrain.all__BMI_Yengo2018.cell_type_results.txt"
# file.ldsc_cts <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.tabula_muris.all__BMI_Yengo2018.cell_type_results.txt"
# file.ldsc_cts <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.campbell_lvl1.all__BMI_Yengo2018.cell_type_results.txt"
# file.ldsc_cts <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.campbell_lvl2.all__BMI_Yengo2018.cell_type_results.txt"
df.ldsc_cts <- load_ldsc_cts_results(file.ldsc_cts, dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")


### Load - MULTI GWAS
dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern=sprintf("%s.(.*).cell_type_results.txt", genomic_annotation_prefix))
filenames <- filenames[!grepl(pattern="__CONDITIONAL__", filenames, perl=T)] # exclude any conditional results
list.dfs <- lapply(file.path(dir.data, filenames), load_ldsc_cts_results, dataset_prefix)
names(list.dfs) <- stringr::str_match(filenames, pattern=sprintf("%s__(.*).cell_type_results.txt", genomic_annotation_prefix))[,2] # ALT: filenames
names(list.dfs)
df.ldsc_cts <- list.dfs %>% bind_rows(.id="gwas")
df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")

### format file [*not* compatible with multi-GWAS plotting]
df.ldsc_cts.export <- df.ldsc_cts %>% select(gwas, p.value, annotation) %>% spread(key=gwas, value=p.value) # format file

# ======================================================================= #
# =============================== METADATA ================================= #
# ======================================================================= #

df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts.export <- df.ldsc_cts.export %>% left_join(df.metadata, by="annotation") # add meta data

# ======================================================================= #
# =============================== EXPORT ================================= #
# ======================================================================= #

### Multi GWAS
df.ldsc_cts.export %>% write_csv(sprintf("out.tmp.190111.%s.sem_mean.multi_gwas.csv", dataset_prefix))

### Single GWAS
# df.ldsc_cts %>% write_csv(sprintf("out.tmp.190111.%s.sem_mean.%s.csv", dataset_prefix, "BMI_UPDATE_Yengo2018"))

# ======================================================================= #
# ================================ PLOT: cell priori [SINGLE-GWAS] ================================= #
# ======================================================================= #


df.plot <- df.ldsc_cts
fdr_threshold <- 0.05/nrow(df.plot)
p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
  geom_point(aes(color=color_by_variable, size=estimate)) + 
  ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
  geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
  labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
  ggtitle(sprintf("%s - %s", "BMI_Yengo2018", dataset_prefix)) +
  theme_classic() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p
file.out <- sprintf("out.tmp.190103.plot.cell_prioritization.%s.%s.sem_mean.color_by_class.pdf", dataset_prefix, "BMI_UPDATE_Yengo2018")
ggsave(p, filename=file.out, width=20, height=8)



# ======================================================================= #
# ================================ PLOT: cell priori [MULTI-GWAS] ================================= #
# ======================================================================= #


list.dfs.split <- split(df.ldsc_cts, df.ldsc_cts$gwas)
for (name_elem in names(list.dfs.split)) {
  print(name_elem)
  df.plot <- list.dfs.split[[name_elem]]
  
  fdr_threshold <- 0.05/nrow(df.plot)
  p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
    geom_point(aes(color=color_by_variable, size=estimate)) + 
    ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
    geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
    labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
    # ggtitle(sprintf("%s - %s", "BMI_Yengo2018", name_elem)) +
    ggtitle(sprintf("%s - %s", dataset_prefix, name_elem)) +
    theme_classic() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  file.out <- sprintf("out.tmp.190111.plot.cell_prioritization.%s.%s.sem_mean.color_by_class.pdf", dataset_prefix, name_elem)
  ggsave(p, filename=file.out, width=20, height=8)
}


# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #



# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #




# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

