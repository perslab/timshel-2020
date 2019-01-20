############### SYNOPSIS ###################
# Postanalysis of LDSC on WGCNA modules

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/"
setwd(wd)


library(tidyverse)


# ======================================================================= #
# ============================ LOAD DATA ============================== #
# ======================================================================= #

dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"

filenames <- list.files(path=dir.data,  pattern="wgcna.maca.*cell_type_results.txt") 
filenames
list.dfs <- lapply(file.path(dir.data, filenames), read_tsv)
names(list.dfs) <- stringr::str_match(filenames, pattern="wgcna.maca.(.*).cell_type_results.txt")[,2] # ALT: filenames
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
# file.multi_geneset <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv" # novo_lira_sema
# file.multi_geneset <- "/projects/jonatan/tmp-mousebrain/tables/mousebrain_Neurons_ClusterName_2_cell_cluster_module_genes.csv" # mousebrain neurons
file.multi_geneset <- "/projects/jonatan/tmp-maca/tables/maca_tissue_cell_type_kME_cell_cluster_module_genes.csv" # maca
df.multi_geneset <- read_csv(file.multi_geneset)
df.multi_geneset <- df.multi_geneset %>% group_by(module) %>% summarise(cell_cluster=unique(cell_cluster), NGENES_MOUSE=n_distinct(ensembl))
# df.multi_geneset <- df.multi_geneset %>% select(cell_cluster,module) %>% distinct() # reduce to unique origin/module_id rows | ALT WORKS

### Join
df <- df %>% left_join(df.multi_geneset, by=c("module_id"="module"))
df <- df %>% rename(module_origin = cell_cluster)

### add more cols
df <- df %>% mutate(module_name=paste(module_origin, module_id, sep="-"))
df$tissue <- stringr::str_split_fixed(df$module_origin,pattern="_",n=2)[,1] # only for MACA
df$cell_type <- stringr::str_split_fixed(df$module_origin,pattern="_",n=2)[,2] # only for MACA

# ======================================================================= #
# ============================ WRITE 'META-DATA FILE' =============================== #
# ======================================================================= #

df.metadata <- df %>% select(gwas, module_id, module_ldsc_pval=P, module_origin,
                             cell_type, tissue, n_genes_module_mouse=NGENES_MOUSE)

### WRITE 'META-DATA FILE'
# df.metadata %>% write_csv("/raid5/projects/timshel/sc-genetics/sc-genetics/out/module_correlation_graph/tabula_muris.module_prioritization.ldsc.metadata.MULTI_TRAIT.csv")

# ======================================================================= #
# ============================ TRAIT SPECIFICITY =============================== #
# ======================================================================= #

df.trait_specificity <- df %>% select(gwas, module_id, module_origin, P) %>% spread(key=gwas, value=P) 
df.trait_specificity %>% write_csv("trait_specificity.maca.csv")
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

df.select <- df.select %>% arrange(tissue) %>% mutate(module_name=factor(module_name, levels=unique(module_name))) # ordering




### facet wrap
p <- ggplot(df.select %>% filter(-log10(P) > 2), aes(x=module_name, y=-log10(P), fill=tissue)) + 
  geom_col() + 
  # facet_wrap(~gwas, ncol = 1) + 
  # guides(fill=FALSE) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
p
# ggsave("out.ldsc.wgcna_maca.barplot_pvals.filter.pdf", w=12, h=10)
# ggsave("out.ldsc.wgcna_maca.barplot_pvals.pdf", w=25, h=10)



# ======================================================================= #
# ======================== Compare LDSC and MAGMA ======================= #
# ======================================================================= #

df.cmp.ldsc <- read_csv("/raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/out.ldsc.wgcna_maca.BMI_Yengo2018.csv")
df.cmp.magma <- read_csv("/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-geneset/out.magma_geneset.wgcna_maca.BMI_Yengo2018.csv")
df.cmp <- full_join(df.cmp.ldsc, df.cmp.magma, by="module_id", suffix=c(".ldsc", ".magma"))
df.cmp
ggplot(df.cmp, aes(x=-log10(P.ldsc), y=-log10(P.magma), color=tissue)) + geom_point(alpha=0.3) + guides(color=FALSE)
with(df.cmp, cor.test(-log10(P.ldsc), -log10(P.magma), method="spearman"))

# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #



# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #
