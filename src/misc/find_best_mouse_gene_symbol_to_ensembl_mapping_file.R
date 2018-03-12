# AIM: identify best mapping file.

library(tidyverse)
library(stringr)

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/misc/"
setwd(wd)

# ======================================================================= #
# ==========================  LOAD MAPPING FILES  ======================== #
# ======================================================================= #


### PT ensembl
file.mapping.ens <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
df.mapping.ens <- read_tsv(file.mapping.ens)

### TP NCBI
file.mapping.ncbi <- "/data/genetic-mapping/ncbi/Mus_musculus.gene_info_symbol2ensembl.gz"
df.mapping.ncbi <- read_tsv(file.mapping.ncbi)
file.ncbi_raw <- "/data/genetic-mapping/ncbi/Mus_musculus.gene_info.gz"
df.ncbi_raw <- read_tsv(file.ncbi_raw)
# SEE also: /projects/tp/tmp-bmi-brain/src/generate_gene_synonym_mapping_mmusculus.R
# README: Files Mus_musculus.gene_info.gz and Homo_sapies.gene_info.gz downloaded from ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/ at March 7, 2018.

### Orthologs
file.ortholog_mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
df.ortholog_mapping <- read_tsv(file.ortholog_mapping)

# ======================================================================= #
# ============================  LOAD GENES  ============================== #
# ======================================================================= #

# save(genes.maca, genes.lira, file="tmp.maca_lira_genes.RData")
load(file="tmp.maca_lira_genes.RData") # genes.maca, genes.lira

### LIRA
# file.RData.cell_atlas <- "/projects/timshel/sc-arc_lira/src/out-generate_cell_atlas/arc_lira_cell_atlas-core_objs.RData"
# load(file.RData.cell_atlas) # seurat_obj, df.cluster_markers
# genes.lira <- seurat_obj@data@Dimnames[[1]]
length(genes.lira) # 20029

### MACA
# file.RData.cell_atlas <- "/projects/timshel/maca/data/figshare/180126-facs/maca.seurat_obj.facs.figshare_180126.RData" # 6.1 GB | seurat_obj
# load(file.RData.cell_atlas) 
# genes.maca <- seurat_obj@data@Dimnames[[1]]
length(genes.maca) # 22977



# =========================================================================================== #
# ================================  Function mapping rate  =================================== #
# ============================================================================================ #


mapping_rate <- function(gene_symbols, df.gene_name_mapping, colname_gene_symbol, colname_ensembl) {
  
  ### for testing
  #gene_symbols=genes.lira
  #df.gene_name_mapping=df.mapping.ens
  #colname_gene_symbol="gene_name_optimal"
  #colname_ensembl="ensembl_gene_id"
  
  df.gene_name_mapping <- as.data.frame(df.gene_name_mapping) # to ensure that we can use data frame operations/indexing that might fail on tibbles
  
  # ==========================  Map: MGI to EnsemblID  ==================== #
  ### Genes not found
  bool.not_found <- !(gene_symbols %in% df.gene_name_mapping[,colname_gene_symbol]) 
  not_mapped_from_MGI_to_ensembl.n_genes <- sum(bool.not_found) # n genes could not be mapped
  not_mapped_from_MGI_to_ensembl.genes <- gene_symbols[bool.not_found] # e.g Fam150a not found, but this is because it is called ALKAL1. Fam150a is listed as a synonym in Emsembl now.
  
  ### Match
  idx.match <- match(gene_symbols, df.gene_name_mapping[,colname_gene_symbol]) # find mapping idx. 'nomatch' is set to NA
  genes.ensembl <- as.character(df.gene_name_mapping[idx.match,colname_ensembl]) # vector with same length as "gene_symbols", but we NA values for genes not found.
  genes.ensembl <- genes.ensembl[!is.na(genes.ensembl)] # removing na values (unmapped genes)
  
  # ========================  Subset on orthologs  ======================== #
  ### Genes without orthologs
  bool.not_found <- !(genes.ensembl %in% df.ortholog_mapping$mmusculus_homolog_ensembl_gene) 
  not_mapped_to_human_ortholog.n_genes <- sum(bool.not_found) # genes with no ortholog
  not_mapped_to_human_ortholog.genes <- genes.ensembl[bool.not_found] # genes with no ortholog
  
  
  # ===============================  PRINT  ============================== #
  
  str1 <- sprintf("Number of genes not mapped from MGI to Ensembl ID: %s out of %s genes (%.2f pct)", 
                  not_mapped_from_MGI_to_ensembl.n_genes, 
                  length(gene_symbols),
                  not_mapped_from_MGI_to_ensembl.n_genes/length(gene_symbols)*100)
  str2 <- sprintf("Number of genes not mapped to human ortholog: %s out of %s genes (%.2f pct)", 
                  not_mapped_to_human_ortholog.n_genes,
                  length(genes.ensembl),
                  not_mapped_to_human_ortholog.n_genes/length(genes.ensembl)*100)
  str3 <- sprintf("Overall mapping rate: %s not mapped out of %s genes (%.2f pct)", 
                  not_mapped_from_MGI_to_ensembl.n_genes+not_mapped_to_human_ortholog.n_genes,
                  length(gene_symbols),
                  (not_mapped_from_MGI_to_ensembl.n_genes+not_mapped_to_human_ortholog.n_genes)/length(gene_symbols)*100)
  print(str1)
  print(str2)
  print(str3)
  
}


# =========================================================================================== #
# ============================  Run mapping rate  ============================== #
# ============================================================================================ #

### Lira - PT ens
mapping_rate(gene_symbols=genes.lira, df.gene_name_mapping=df.mapping.ens, colname_gene_symbol="gene_name_optimal", colname_ensembl="ensembl_gene_id")
### Lira - TP NCBI
mapping_rate(gene_symbols=genes.lira, df.gene_name_mapping=df.mapping.ncbi, colname_gene_symbol="symbol", colname_ensembl="ensembl")

### MACA - PT ens
mapping_rate(gene_symbols=genes.maca, df.gene_name_mapping=df.mapping.ens, colname_gene_symbol="gene_name_optimal", colname_ensembl="ensembl_gene_id")
### MACA - TP NCBI
mapping_rate(gene_symbols=genes.maca, df.gene_name_mapping=df.mapping.ncbi, colname_gene_symbol="symbol", colname_ensembl="ensembl")


    
# =========================================================================================== #
# ============================  Check for uniqueness of symbols ============================== #
# ============================================================================================ #



### TP NCBI file (raw)
nrow(df.ncbi_raw) # 68844
n_distinct(df.ncbi_raw$Symbol_from_nomenclature_authority) # 65649
n_distinct(df.ncbi_raw$Symbol) # 68506
# --> symbols are not unique, which is likely an issue


### TP NCBI file
nrow(df.mapping.ncbi) # 86602
n_distinct(df.mapping.ncbi$ensembl) # 24587
# --> few ensembl ids is not good


### PT Ens file
nrow(df.mapping.ens) # 53443
n_distinct(df.mapping.ens$gene_name_optimal) # 53443
n_distinct(df.mapping.ens$ensembl_gene_id) # 52438
# --> almost all entries are unique. All entries would be unique if not for the "optimal mapping feature".

### ...
# PT fil indeholder 52438 EnsemblIDs 
# TP indeholder 24587 EnsemblIDs. 

# =========================================================================================== #
# ============================  Check for uniqueness of synonyms ============================== #
# ============================================================================================ #

### Make a tidy data frame (split synonyms into multiple rows)
# REF: https://stackoverflow.com/questions/44401023/splitting-multiple-values-in-one-column-into-multiple-rows-r
df.ncbi_raw.s1 <- df.ncbi_raw %>% 
  select(Symbol, Synonyms, dbXrefs) %>%
  separate_rows(Synonyms,sep="\\|")
  
df.ncbi_raw.s2 <- df.ncbi_raw.s1 %>%
  separate_rows(dbXrefs,sep="\\|") %>% 
  filter(grepl("Ensembl", dbXrefs)) %>%
  mutate(ensembl_id = sapply(str_split(dbXrefs,":"), '[[', 2)) # extracting "ENSMUSG00000030359" from "Ensembl:ENSMUSG00000030359"

df.ncbi_raw.tidy <- as.data.frame(df.ncbi_raw.s2) # copy

### Check that synonyms are not unique
multimapping_gene_symbol <- c()
multimapping_synonyms <- c()
idx_progress <- 0
for (ensembl_id in unique(df.ncbi_raw.tidy$ensembl_id)) {
  idx_progress <- idx_progress + 1 
  print(sprintf("Processing %s/%s", idx_progress, n_distinct(df.ncbi_raw.tidy$ensembl_id)))
  df.current_ensembl_id <- df.ncbi_raw.tidy[df.ncbi_raw.tidy$ensembl_id==ensembl_id,]
  df.other_ensembl_id <- df.ncbi_raw.tidy[df.ncbi_raw.tidy$ensembl_id!=ensembl_id,]
  current_gene_symbols <- unique(df.current_ensembl_id$Symbols)
  for (gene_symbol in current_gene_symbols) {
    if (any((gene_symbol %in% df.other_ensembl_id$Symbol))) {
      multimapping_gene_symbol <- unique(c(multimapping_gene_symbol, gene_symbol))
      print(sprintf("N multimapping_gene_symbol = %s", length(multimapping_gene_symbol)))
    }
  }
  for (gene_synonym in df.current_ensembl_id$Synonyms) {
    if (any((gene_synonym %in% df.other_ensembl_id$Synonyms))) {
      multimapping_synonyms <- unique(c(multimapping_synonyms, gene_synonym))
      print(sprintf("N multimapping_synonyms = %s", length(multimapping_synonyms)))
    }
  }
}
length(multimapping_gene_symbol) # 
length(multimapping_synonyms) # many multimapping
