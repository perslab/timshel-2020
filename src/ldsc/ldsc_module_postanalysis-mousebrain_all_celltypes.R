############### SYNOPSIS ###################
# Postanalysis of LDSC on WGCNA modules from Mousebrain ALL CELL-TYPES

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/ldsc/"
setwd(wd)


library(tidyverse)

# ======================================================================= #
# ================= Prepare combined .multi_geneset file =============== #
# ======================================================================= #

### ONLY NEED TO DO THIS ONCE (hopefully)

# files.relative <- c("mousebrain_Astrocytes/log.mousebrain_Astrocytes.multi_geneset.txt",
#                     "mousebrain_Ependymal/log.mousebrain_Ependymal.multi_geneset.txt",
#                     "mousebrain_Immune/log.mousebrain_Immune.multi_geneset.txt",
#                     "mousebrain_Neurons/log.Neurons_ClusterName.multi_geneset.txt", # OBS: different format for naming
#                     "mousebrain_Oligos/log.mousebrain_Oligos.multi_geneset.txt",
#                     "mousebrain_PeripheralGlia/log.mousebrain_PeripheralGlia.multi_geneset.txt",
#                     "mousebrain_Vascular/log.mousebrain_Vascular.multi_geneset.txt")
# filenames <- file.path("/scratch/sc-ldsc", files.relative) 
# list.dfs <- lapply(filenames, read_tsv)
# names(list.dfs) <- stringr::str_match(filenames, pattern="/scratch/sc-ldsc/mousebrain_(.*)/log.*.txt")[,2]
# df <- list.dfs %>% bind_rows(.id="wgcna_run")
# names(df)
# 
# ### add new unique identifier
# df <- df %>% mutate(module_color_id = annotation) # rename
# df <- df %>% mutate(annotation = paste0(wgcna_run, ".", module_color_id)) # ---> e.g. Astrocytes.blue1
# 
# # df %>% write_tsv("log.mousebrain_ALL_CLASSES.multi_geneset.txt")

# ======================================================================= #
# ============================ LOAD DATA ============================== #
# ======================================================================= #

dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern="wgcna.mousebrain_(.*).cell_type_results.txt") 
filenames
list.dfs <- lapply(file.path(dir.data, filenames), read_tsv)
names(list.dfs) <- stringr::str_match(filenames, pattern="wgcna.mousebrain_(.*)\\.cell_type_results.txt")[,2] # ALT: filenames
names(list.dfs)
df.bind <- list.dfs %>% bind_rows(.id="origin_and_gwas")

df.ldsc_cts <- df.bind %>% separate(col="origin_and_gwas",into=c("wgcna_run", "gwas"), sep="\\.")
# ======================================================================= #
# ============================ PROCESS =============================== #
# ======================================================================= #

df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=stringr::str_replace_all(Name, pattern="mousebrain_", "")) # SPECIFIC
df.ldsc_cts <- df.ldsc_cts %>% rename(BETA=Coefficient, SE=Coefficient_std_error, P=Coefficient_P_value)
df.ldsc_cts <- df.ldsc_cts %>% group_by(gwas) %>% mutate(P_adj = P*n()) %>% ungroup()


# ======================================================================= #
# ============================ MAP MODULE ORIGIN  =============================== #
# ======================================================================= #
### GENESET
file.module_geneset <- "/scratch/sc-ldsc/log.mousebrain_ALL_CLASSES.multi_geneset.txt"
df.module_geneset <- read_tsv(file.module_geneset)
df.module_geneset <- df.module_geneset %>% rename(ensembl_gene_id=gene, module_id=annotation, module_origin=cell_cluster)
df.module_metadata <- df.module_geneset %>% select(module_origin, module_id) %>% distinct()

df <- df.ldsc_cts %>% select(gwas, module_id, module_ldsc_pval=P) %>% left_join(df.module_metadata, by="module_id")

# ======================================================================= #
# ============================ ADD CELL META DATA  =============================== #
# ======================================================================= #

file.module_origin_metadata <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-agg_L5.metadata.csv"
df.module_origin_metadata <- read_csv(file.module_origin_metadata)
df.module_origin_metadata <- df.module_origin_metadata %>% select(module_origin=ClusterName,
       Developmental_compartment,
       Probable_location,
       Region,
       TaxonomyRank1,
       TaxonomyRank2,
       TaxonomyRank3,
       TaxonomyRank4)
       
df.metadata <- df %>% left_join(df.module_origin_metadata, by="module_origin")

### WRITE 'META-DATA FILE'
# df.metadata %>% write_csv("mousebrain.module_prioritization.ldsc.metadata.MULTI_TRAIT.csv")


# ======================================================================= #
# ============================ TRAIT SPECIFICITY =============================== #
# ======================================================================= #

df.trait_specificity <- df %>% spread(key=gwas, value=module_ldsc_pval) 
df.trait_specificity %>% write_csv("trait_specificity.mousebrain.csv")

# df.trait_specificity %>% write_csv("trait_specificity.hypothalamus_mette_thesis.csv")

df.x <- df.trait_specificity %>% select(-module_id, -module_origin)
df.cor <- cor(-log10(df.x))












# ======================================================================= #
# ============================ FILTER =============================== #
# ======================================================================= #

df.bmi <- df %>% filter(gwas=='BMI_Yengo2018') %>% arrange(P)
# df.bmi %>% write_csv("out.ldsc.wgcna_maca.BMI_Yengo2018.csv")

df.select <- df.bmi # SWITCH

# ======================================================================= #
# ============================ PLOT =============================== #
# ======================================================================= #

df.select <- df.select %>% arrange(module_origin) %>% mutate(module_name=factor(module_name, levels=unique(module_name))) # ordering

### facet wrap
p <- ggplot(df.select, aes(x=fct_reorder(module_name, -log10(P)), y=-log10(P), fill=module_origin)) + # %>% filter(-log10(P) > 2)
  geom_col() + 
  # facet_wrap(~gwas, ncol = 1) + 
  guides(fill=FALSE) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
p
# ggsave("out.ldsc.wgcna_maca.barplot_pvals.filter.pdf", w=12, h=10)
# ggsave("out.ldsc.wgcna_maca.barplot_pvals.pdf", w=25, h=10)


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #



# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #
