############### SYNOPSIS ###################
# Postanalysis of LDSC on WGCNA modules from Mousebrain NEURONS.

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
# ============================ LOAD DATA ============================== #
# ======================================================================= #

dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"

filenames <- list.files(path=dir.data,  pattern="wgcna.mousebrain.*cell_type_results.txt") 
# filenames <- list.files(path=dir.data,  pattern="wgcna.hypothalamus_mette_thesis.*cell_type_results.txt") # METTE
filenames
list.dfs <- lapply(file.path(dir.data, filenames), read_tsv)
names(list.dfs) <- stringr::str_match(filenames, pattern="wgcna.mousebrain.(.*).cell_type_results.txt")[,2] # ALT: filenames
# names(list.dfs) <- stringr::str_match(filenames, pattern="wgcna.hypothalamus_mette_thesis.(.*).cell_type_results.txt")[,2] # ALT: filenames
names(list.dfs)
df <- list.dfs %>% bind_rows(.id="gwas")


# ======================================================================= #
# ============================ PROCESS =============================== #
# ======================================================================= #

### process
df$module_id <- stringr::str_split_fixed(df$Name,pattern="\\.",n=2)[,2]
### ALT: this also works, but the 'module_origin' is not informative with the current CTS files.
# df <- df %>% separate(Name, into=c("module_origin", "module_id"), sep="\\.(?=[^.]*$)", extra="merge", remove=F) # sep is regex and matches LAST DOT in string. *NOT sure about the escaping*
# ^ REGEX REF DOT: https://stackoverflow.com/a/19976308/6639640
# ^ the above separate() will split on the LAST DOT (And split in only two parts)
# ^ REGEX REF: https://stackoverflow.com/a/11134049/6639640: Matches a literal hypthen and then asserts ((?=) is positive lookahead) that no other character up to the end of the string is a hypthen 

### rename to match MAGMA
df <- df %>% rename(BETA=Coefficient, SE=Coefficient_std_error, P=Coefficient_P_value)

### Calc adjust P-val
df <- df %>% group_by(gwas) %>% mutate(P_adj = P*n())


# ======================================================================= #
# ============================ MAP module origin =============================== #
# ======================================================================= #

### Read geneset file
file.multi_geneset <- "/projects/jonatan/tmp-mousebrain/tables/mousebrain_Neurons_ClusterName_2_cell_cluster_module_genes.csv" # mousebrain neurons
# file.multi_geneset <- "/projects/jonatan/tmp-maca/tables/maca_tissue_cell_type_kME_cell_cluster_module_genes.csv" # maca
# file.multi_geneset <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_lists/mludwig_thesis_hypothalamus_wgcna_modules.csv" # mette hypo
df.multi_geneset <- read_csv(file.multi_geneset)
df.multi_geneset <- df.multi_geneset %>% group_by(module) %>% summarise(cell_cluster=unique(cell_cluster), NGENES_MOUSE=n_distinct(ensembl))
# df.multi_geneset <- df.multi_geneset %>% select(cell_cluster,module) %>% distinct() # reduce to unique origin/module_id rows | ALT WORKS

### Join
df <- df %>% left_join(df.multi_geneset, by=c("module_id"="module"))
df <- df %>% rename(module_origin = cell_cluster)

### add more cols
df <- df %>% mutate(module_name=paste(module_origin, module_id, sep="-"))
# df$tissue <- stringr::str_split_fixed(df$module_origin,pattern="_",n=2)[,1] # only for MACA
# df$cell_type <- stringr::str_split_fixed(df$module_origin,pattern="_",n=2)[,2] # only for MACA


# ======================================================================= #
# ============================ TRAIT SPECIFICITY =============================== #
# ======================================================================= #

df.trait_specificity <- df %>% select(gwas, module_id, module_origin, P) %>% spread(key=gwas, value=P) 
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
