### SYNOPSIS: read SEM data; map ensembl to entrez; write output files.

library(tidyverse)

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/"
setwd(wd)

# ======================================================================= #
# =============================== FUNCTIONS ================================ #
# ======================================================================= #

### Map:  ENSEMBL --> ENTREZ
map_ensembl2entrez <- function(df) {
  ### INPUT df: a tibble/data.frame with the column 'gene' with ensembl gene ids.
  ### OUTOUT df: a tibble with genes mapped to entrez ids. Only mapped genes are retained.
  
  ### TODO: make 'gene' column a variable, so the name does not have to be 'gene'
  ### TODO: check the entrez mapped genes are UNIQUE
  
  file.mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_id_mapping.hsapiens.ensembl_entrez.txt.gz"
  df.mapping <- read_tsv(file.mapping) # col1=ensembl_gene_id, col2=entrezgene
  
  genes_mapped <- df.mapping$entrezgene[match(df$gene, df.mapping$ensembl_gene_id)]
  bool_dups <- duplicated(na.omit(genes_mapped)) # excluding NA when counting duplicated. duplicated(c(NA,NA,NA)) returns FALSE  TRUE  TRUE.
  print(sprintf("Number of genes without unique mapping to entrez (genes with duplicated entrez ids after mapping): %s",sum(bool_dups)))
  print(sprintf("Number of genes mapped: %s",sum(!is.na(genes_mapped))))
  print(sprintf("Number of genes not mapped: %s",sum(is.na(genes_mapped)))) # number of not mapped genes
  df.clean <- df %>% 
    mutate(gene=genes_mapped) %>%
    filter(!is.na(gene)) %>%
    filter(!duplicated(gene)) # keep only one of the duplicated pair (if any)
  return(df.clean)
}

# ============================================================================ #
# ======================= READ SEM DATA ================ #
# ============================================================================ #

### read SEM expression data
DIR.sem_data <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain"
filenames <- list.files(path=DIR.sem_data,  pattern="*.hsapiens_orthologs.csv.gz") # has *HUMAN* ENSEMBL IDs
list.df_sem <- lapply(file.path(DIR.sem_data, filenames), read_csv)
filenames_shorten <- stringr::str_match(filenames, "celltype_expr\\.(.*)\\.hsapiens_orthologs.csv.gz")[,2] # e.g. "celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz" --> "specificity_quantiles"
names(list.df_sem) <- filenames_shorten
names(list.df_sem)

### testing
# df.x <- map_ensembl2entrez(list.df_sem[[1]])

### map to entrez
list.df_sem.mapped <- lapply(list.df_sem, map_ensembl2entrez)
# [1] "Number of genes without unique mapping to entrez (genes with duplicated entrez ids after mapping): 1"
# [1] "Number of genes mapped: 14813"
# [1] "Number of genes not mapped: 132"

### write output files
for (name in names(list.df_sem.mapped)) {
  print(name)
  file.out <- sprintf("mousebrain.sem.%s.entrez.txt", name)
  write_tsv(list.df_sem.mapped[[name]], file.out)
}


