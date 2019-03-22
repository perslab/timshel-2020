############### SYNOPSIS ###################
# Add mouse ortholog genes to MAGMA .genes.out file.

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================= #
# =============================== FUNCTIONS ================================ #
# ======================================================================= #

### Map:  HUMAN --> MOUSE
human_to_mouse_ortholog_gene_expression_mapping <- function(df, colname_geneids_from="ensembl_gene_id", colname_geneids_to="mmusculus_homolog_ensembl_gene") {
  ### INPUT df: a tibble/data.frame with the column 'colname_geneids_from' with human ensembl gene ids.
  ### OUTOUT df: a tibble with MOUSE ensembl gene ids added to the column 'colname_geneids_to'. Genes that did not map have NA values.
  
  file.mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
  df.mapping <- suppressMessages(read_delim(file.mapping, delim="\t"))
  # ensembl_gene_id chromosome_name start_position  end_position    mmusculus_homolog_ensembl_gene  mmusculus_homolog_orthology_confidence
  # ENSG00000138593 15      49280673        49338760        ENSMUSG00000035093      1
  # ENSG00000166351 21      14982498        15013906        ENSMUSG00000095294      0
  df <- as.data.frame(df) # convert to df to ensure the the below operations work.
  genes_mapped <- df.mapping$mmusculus_homolog_ensembl_gene[match(df[,colname_geneids_from], df.mapping$ensembl_gene_id)]
  bool_dups <- duplicated(na.omit(genes_mapped)) # excluding NA when counting duplicated. duplicated(c(NA,NA,NA)) returns FALSE  TRUE  TRUE.
  print(sprintf("Number of genes with a NON-unique mapping (genes with duplicated gene IDs after mapping): %s",sum(bool_dups)))
  print(sprintf("Number of genes mapped: %s",sum(!is.na(genes_mapped))))
  print(sprintf("Number of genes not mapped: %s",sum(is.na(genes_mapped)))) # number of not mapped genes
  df <- df %>% mutate(!!rlang::sym(colname_geneids_to):=genes_mapped)
  # filter(!is.na(gene)) %>% # remove all rows without mapping
  # filter(!duplicated(gene)) # keep only one of the duplicated pair (if any)
  return(df)
}
  


# ======================================================================= #
# ======================== READ MAGMA .genes.out ======================== #
# ======================================================================= #
# file.magma <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/BMI_Yengo2018.correct_all.resid.genes.out"
file.magma <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/BMI_Yengo2018.resid.correct_all.gsa.genes.out"
df.magma <- read_table(file.magma, comment = "#")

### add ensembl ids
df.magma <- add_ensembl_ids_from_entrez(df.magma, colname_geneids_from="GENE", colname_geneids_to="ensembl_gene_id")
# [1] "Number of genes mapped: 17467"
# [1] "Number of genes not mapped: 158"
# [1] "Number of genes with a NON-unique mapping (genes with duplicated ensembl gene IDs after mapping): 33"
# [1] "Total mapping stats: 191 genes have no mapping (not mapped + duplicates) out of 17625 input genes."
# [1] "Total genes mapped (non NA genes): 17434"
# [1] "Returning tibble with the column 'ensembl_gene_id' added where all gene identifiers unique. Unmapped genes have NA values"

# ======================================================================= #
# =========================== Map from human to mouse =================== #
# ======================================================================= #

### add mouse orthologs
df.magma <- human_to_mouse_ortholog_gene_expression_mapping(df.magma, colname_geneids_from="ensembl_gene_id", colname_geneids_to="mmusculus_homolog_ensembl_gene")
# [1] "Number of genes with a NON-unique mapping (genes with duplicated gene IDs after mapping): 0"
# [1] "Number of genes mapped: 14782"
# [1] "Number of genes not mapped: 2843"


sum(is.na(df.magma$mmusculus_homolog_ensembl_gene))


# ======================================================================= #
# ============================= Write out file ================================ #
# ======================================================================= #

# df.magma %>% write_tsv("BMI_Yengo2018.resid.correct_all.mouse_homologs.gsa.genes.out")




# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #







