############### SYNOPSIS ###################
### Investigate MAGMA P-value bias


library(tidyverse)
library(ggpubr) # http://www.sthda.com/english/rpkgs/ggpubr/index.html


# ========================================================================== #
# ======================= SETUP ======================= #
# ========================================================================== #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/magma/"
setwd(wd)

# ========================================================================== #
# ======================= READ MAGMA GWAS PVALUES FILES ======================= #
# ========================================================================== #


### GET MAPPING: ENTREZ --> ENSEMBL
ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
df.mapping <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'entrezgene'), mart = ensembl) # data frame
head(df.mapping)


### READ MAGMA PVALUES
DIR.gwas <- "/scratch/tmp-magma_gwas/"
filenames <- list.files(DIR.gwas, pattern="*genes.out")
list.magma <- lapply(file.path(DIR.gwas, filenames), read_table) # read_table() and read_table2() are designed to read the type of textual data where each column is separated by one (or more) columns of space.
names(list.magma) <- stringr::str_match(filenames, pattern="(.*).txt*")[,2] # e.g. "AN_PGC_Duncan2017.txt.10UP.1.5DOWN.genes.out" --> "AN_PGC_Duncan2017"
names(list.magma)

### Map:  ENTREZ --> ENSEMBL
map_entrez2ensembl <- function(df.magma) {
  genes_mapped <- df.mapping$ensembl_gene_id[match(df.magma$GENE, df.mapping$entrezgene)]
  print(sum(is.na(genes_mapped))) # number of not mapped genes
  df.magma.clean <- df.magma %>% 
    mutate(ensembl_gene_id=genes_mapped) %>%
    filter(!is.na(ensembl_gene_id))
  return(df.magma.clean)
}

list.magma <- lapply(list.magma, map_entrez2ensembl) # run function

### combine list of dfs to a single df
df.magma <- bind_rows(list.magma, .id="gwas")

### compute length
df.magma <- df.magma %>% mutate(LENGTH=STOP-START)


# ========================================================================== #
# ======================= INVESTIGATE BIAS ======================= #
# ========================================================================== #

# gene_length vs ZSTAT
# NSNPS vs gene_length
# Spearmann correlation: rnorm() vs ZSTAT


df.res <- df.magma %>% 
  group_by(gwas) %>%
  summarise(cor_p_length_vs_zstat = cor.test(LENGTH, ZSTAT)$p.value, 
            cor_p_length_vs_nsnps = cor.test(LENGTH, NSNPS)$p.value,
            cor_p_zstat_vs_rnorm = cor.test(ZSTAT, rnorm(n()), method="spearman")$p.value)


### Plot: LENGTH VS ZSTAT
p <- ggplot(df.magma, aes(x=LENGTH, y=ZSTAT)) + geom_point(alpha=0.2) + geom_smooth(method = "lm", se = FALSE) + facet_wrap(~gwas)
p <- p + ggpubr::stat_cor(method = "spearman")
# p
# ggsave("plot.magma_length_bias.pdf", w=30, h=30)



