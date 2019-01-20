############### SYNOPSIS ###################
# Postanalysis of MAGMA gene set test (on WGCNA modules)

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-geneset/"
setwd(wd)


library(tidyverse)


### SNIPPET *.sets.out
# # MEAN_SAMPLE_SIZE = 137968
# # TOTAL_GENES = 17728
# # CONDITIONED_INTERNAL = genesize, log_genesize, genedensity, log_genedensity, samplesize, log_samplesize, inverse_mac, log_inverse_mac
# SET              NGENES       BETA   BETA_STD         SE            P
# darkkhaki            27    -0.0249   -0.00097      0.157      0.56297
# darkmagenta          61    -0.0485   -0.00284      0.105      0.67834
# maroon              110     0.0526    0.00413     0.0811       0.2584
# mediumturquoise      84     0.0181    0.00124        0.1      0.42856
# skyblue             685    -0.0249   -0.00479     0.0321      0.78093


# ======================================================================= #
# ============================ LOAD DATA ============================== #
# ======================================================================= #

dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-geneset/out.maca/"

filenames <- list.files(path=dir.data,  pattern=".*.sets.out") 
filenames
list.dfs <- lapply(file.path(dir.data, filenames), read_table, comment="#")
names(list.dfs) <- stringr::str_match(filenames, pattern="magma_geneset_(.*).sets.out")[,2] # ALT: filenames
names(list.dfs)
df <- list.dfs %>% bind_rows(.id="gwas")

### process
if ("FULL_NAME" %in% colnames(df)) { # if MAGAM abreviated the gene set
  df <- df %>% mutate(SET = FULL_NAME)
  df <- df %>% select(-FULL_NAME)
}
df <- df %>% mutate(SET = stringr::str_replace_all(SET, ":WS:", " ")) # substitute "placeholder whitespace" with whitespace
df <- df %>% separate(SET, into=c("module_origin", "module_id"), sep="-(?=[^-]*$)", extra="merge", remove=F) # sep is regex and matches LAST hyphen in string.
# ^ the above separate() will split on the LAST hyphen. (And split in only two parts)
# ^ REGEX REF: https://stackoverflow.com/a/11134049/6639640: Matches a literal hypthen and then asserts ((?=) is positive lookahead) that no other character up to the end of the string is a hypthen 
df

### ***only for MACA***
df$tissue <- stringr::str_split_fixed(df$module_origin,pattern="_",n=2)[,1] # only for MACA
df$cell_type <- stringr::str_split_fixed(df$module_origin,pattern="_",n=2)[,2] # only for MACA

# ======================================================================= #
# ============================ FILTER =============================== #
# ======================================================================= #

df.bmi <- df %>% filter(gwas=='BMI_Yengo2018') %>% arrange(P)
# df.bmi %>% write_csv("out.magma_geneset.wgcna_maca.BMI_Yengo2018.csv")

df.bmi.hbcho3 <- df %>% filter(gwas=='BMI_Yengo2018', module_origin=='HBCHO3')

# ======================================================================= #
# ============================ PLOT: P vs NGENES =============================== #
# ======================================================================= #

ggplot(df, aes(x=NGENES, y=-log10(P), color=module_origin)) + geom_point(alpha=0.3) + guides(color=FALSE)
with(df, cor.test(NGENES, -log10(P), method="spearman"))
# S = 2.8613e+12, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.1157848 ---> more genes in geneset ==> higher significance


# ======================================================================= #
# ============================ PLOT =============================== #
# ======================================================================= #

### BMI MACA PLOT
df.bmi <- df.bmi %>% arrange(tissue) %>% mutate(SET=factor(SET, levels=unique(SET))) # ordering
p <- ggplot(df.bmi, aes(x=SET, y=-log10(P), fill=tissue)) + 
  geom_col() + 
  # guides(fill=FALSE) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
# ggsave("out.magma_geneset.wgcna_maca.barplot_pvals.pdf", w=25, h=10)


### facet wrap
ggplot(df, aes(x=SET, y=-log10(P), fill=module_origin)) + geom_col() + facet_wrap(~gwas, ncol = 1) + theme(axis.text.x=element_text(angle = 45, hjust = 1))
# ggsave("plot.magma_module_geneset_test.barplot_multigwas.pdf", w=8, h=8)

### dogde (ok, but not great)
# ggplot(df, aes(x=SET, y=-log10(P), fill=gwas)) + geom_col(position=position_dodge()) + theme(axis.text.x=element_text(angle = 45, hjust = 1))




# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #
