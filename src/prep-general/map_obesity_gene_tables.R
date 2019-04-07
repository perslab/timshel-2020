

library(tidyverse)


# ======================= Function ======================= #



add_ensembl_gene_id_from_gene_symbol <- function(df, colname_genes="gene") {
  ### DESCRIPTION: Maps human gene symbols to ensembl gene ids (adds column to data frame).
  
  # ======================================================================= #
  # ================================ INIT  ================================ #
  # ======================================================================= #
  df <- as.data.frame(df) # make sure we don't have a tibble. The below operations (might) fail on tibbles
  
  # ======================================================================= #
  # ==============================  READ DATA  ============================== #
  # ======================================================================= #
  print("Reading gene mapping files...")
  
  file.gene_name_mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/GRCh38.ens_v90.gene_name_version2ensembl.txt.gz"
  df.gene_name_mapping <- suppressMessages(read_delim(file.gene_name_mapping, delim="\t"))
  
  # ======================================================================= #
  # ==========================  Map: Gene symbol to EnsemblID  ==================== #
  # ======================================================================= #
  print("Mapping from Gene symbol to EnsemblID...")
  
  ### Converting Gene symbols to lowercase to avoid issues with case sensitivity
  print("Converting Gene symbols to lowercase to avoid issues with case sensitivity")
  df.gene_name_mapping$gene_name_optimal <- stringr::str_to_lower(df.gene_name_mapping$gene_name_optimal)
  df[,colname_genes] <- stringr::str_to_lower(df[,colname_genes])
  
  ### Genes not found
  bool.not_found <- !(df[,colname_genes] %in% df.gene_name_mapping$gene_name_optimal) 
  not_mapped_from_symbol_to_ensembl.n_genes <- sum(bool.not_found) # n genes could not be mapped
  not_mapped_from_symbol_to_ensembl.genes <- df[,colname_genes][bool.not_found] # e.g Fam150a not found, but this is because it is called ALKAL1. Fam150a is listed as a synonym in Emsembl now.
  
  ### Match
  idx.match <- match(df[,colname_genes], df.gene_name_mapping$gene_name_optimal) # find mapping idx. 'nomatch' is set to NA
  genes.ensembl <- df.gene_name_mapping$ensembl_gene_id[idx.match] # vector with same length as "df[,colname_genes]", but we NA values for genes not found.
  
  ### Set results
  df.ens <- df
  df.ens$ensembl_gene_id <- genes.ensembl
  df.ens[,colname_genes] <- stringr::str_to_title(df.ens[,colname_genes]) # convert to title case
  
  # ======================================================================= #
  # ===============================  FINISH  ============================== #
  # ======================================================================= #
  str1 <- sprintf("Number of genes not mapped from symbol to Ensembl ID: %s out of %s genes (%.2f pct)", 
                  not_mapped_from_symbol_to_ensembl.n_genes, 
                  length(df[,colname_genes]),
                  not_mapped_from_symbol_to_ensembl.n_genes/length(df[,colname_genes])*100)
  str3 <- sprintf("Output dimension of final data table: n_genes=%s x n_features=%s", 
                  nrow(df.ens),
                  ncol(df.ens))
  print(str1)
  print(str3)
  print("List of genes not mapped from Gene symbol to Ensembl ID:")
  print(not_mapped_from_symbol_to_ensembl.genes)
  
  return(df.ens)
  
}


mouse_mgi_to_human_ortholog_mapping <- function(df, colname_genes_mgi="gene_mgi") {
  ### INPUT
  # df                 data frame to map
  #                         columns: must contain a column with columnname 'colname_mgi' with *MGI symbols*. Other columns allowed and will be returned in output dataframe.
  # file.out.map_stats:     output file path for the mapping stats (by default it just prints to the screen)
  ### OUTPUT
  # df.ens.human:      a data frame where mouse genes have been mapped to human genes (Ensembl)
  
  # ======================================================================= #
  # ================================ INIT  ================================ #
  # ======================================================================= #
  ### for debugging
  # df <- df.genes_mouse_obesity
  # colname_genes_mgi <- "gene"
  
  df <- as.data.frame(df) # make sure we don't have a tibble. The below operations (might) fail on tibbles
  
  # ======================================================================= #
  # ==============================  READ DATA  ============================== #
  # ======================================================================= #
  print("Reading gene mapping files...")
  
  # file.gene_name_mapping <- "../data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
  file.gene_name_mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
  df.gene_name_mapping <- suppressMessages(read_delim(file.gene_name_mapping, delim="\t"))
  
  # file.ortholog_mapping <- "../data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
  file.ortholog_mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
  df.ortholog_mapping <- suppressMessages(read_delim(file.ortholog_mapping, delim="\t"))
  
  # ======================================================================= #
  # ==========================  Map: MGI to EnsemblID  ==================== #
  # ======================================================================= #
  print("Mapping from MGI to EnsemblID...")
  
  ### Converting MGI symbols to lowercase to avoid issues with case sensitivity
  print("Converting MGI symbols to lowercase to avoid issues with case sensitivity")
  df.gene_name_mapping$gene_name_optimal <- stringr::str_to_lower(df.gene_name_mapping$gene_name_optimal)
  df[,colname_genes_mgi] <- stringr::str_to_lower(df[,colname_genes_mgi])
  
  ### Genes not found
  bool.not_found <- !(df[,colname_genes_mgi] %in% df.gene_name_mapping$gene_name_optimal) 
  not_mapped_from_MGI_to_ensembl.n_genes <- sum(bool.not_found) # n genes could not be mapped
  not_mapped_from_MGI_to_ensembl.genes <- df[,colname_genes_mgi][bool.not_found] # e.g Fam150a not found, but this is because it is called ALKAL1. Fam150a is listed as a synonym in Emsembl now.
  
  ### Match
  idx.match <- match(df[,colname_genes_mgi], df.gene_name_mapping$gene_name_optimal) # find mapping idx. 'nomatch' is set to NA
  genes.ensembl <- df.gene_name_mapping$ensembl_gene_id[idx.match] # vector with same length as "df[,colname_genes_mgi]", but we NA values for genes not found.
  
  ### Set results
  df.ens <- df
  df.ens$ensembl_gene_id_mouse <- genes.ensembl
  df.ens[,colname_genes_mgi] <- stringr::str_to_title(df.ens[,colname_genes_mgi]) # convert to title case
  
  ### Filter out genes not mapped, and replace gene names with IDs.
  # df.ens <- df[!is.na(genes.ensembl),] # subset to mapped genes
  # rownames(df.ens) <- genes.ensembl[!is.na(genes.ensembl)] # replace gene names with IDs
  
  # ======================================================================= #
  # ========================  Map orthologs  ======================== #
  # ======================================================================= #
  print("Mapping to orthologs...")
  
  ### Find mouse genes with human orthologs
  idx.match <- match(df.ens$ensembl_gene_id_mouse, df.ortholog_mapping$mmusculus_homolog_ensembl_gene)
  genes.human.orthologs <- df.ortholog_mapping$ensembl_gene_id[idx.match] # NA values for genes with no ortholog.
  
  ### Genes without orthologs
  bool.not_found <- !(df.ens$ensembl_gene_id_mouse %in% df.ortholog_mapping$mmusculus_homolog_ensembl_gene) 
  not_mapped_to_human_ortholog.n_genes <- sum(bool.not_found) # genes with no ortholog
  not_mapped_to_human_ortholog.genes <- df.ens$ensembl_gene_id_mouse[bool.not_found] # genes with no ortholog
  
  df.ens.human <- df.ens
  df.ens.human$ensembl_gene_id_human <- genes.human.orthologs
  ### Subset genes
  # df.ens.human <- df.ens[!is.na(genes.human.orthologs),] # subset to ortholog genes
  # rownames(df.ens.human) <- genes.human.orthologs[!is.na(genes.human.orthologs)] # replace mouse ensembl ID with human
  
  # ======================================================================= #
  # ===============================  FINISH  ============================== #
  # ======================================================================= #
  str1 <- sprintf("Number of genes not mapped from MGI to Ensembl ID: %s out of %s genes (%.2f pct)", 
                  not_mapped_from_MGI_to_ensembl.n_genes, 
                  length(df[,colname_genes_mgi]),
                  not_mapped_from_MGI_to_ensembl.n_genes/length(df[,colname_genes_mgi])*100)
  str2 <- sprintf("Number of genes not mapped to human ortholog: %s out of %s genes (%.2f pct)", 
                  not_mapped_to_human_ortholog.n_genes,
                  length(df.ens$ensembl_gene_id_mouse),
                  not_mapped_to_human_ortholog.n_genes/length(df.ens$ensembl_gene_id_mouse)*100)
  str3 <- sprintf("Output dimension of final data table (ortholog mapped): n_genes=%s x n_features=%s", 
                  nrow(df.ens.human),
                  ncol(df.ens.human))
  print(str1)
  print(str2)
  print(str3)
  print("List of genes not mapped from MGI to Ensembl ID:")
  print(not_mapped_from_MGI_to_ensembl.genes)
  print("List of genes not mapped to human ortholog:")
  print(not_mapped_to_human_ortholog.genes)
  
  return(df.ens.human)
  
}  


# ======================= file.genes_mouse_obesity ======================= #

file.genes_mouse_obesity <- "/projects/timshel/sc-genetics/sc-genetics/data/genes_obesity/original_tables/yazdi2015_table1_mouse_obesity_genes.txt"
df.genes_mouse_obesity <- read_tsv(file.genes_mouse_obesity) 
df.genes_mouse_obesity <- df.genes_mouse_obesity %>% filter(!duplicated(gene)) # remove dups
### HGNChelper does not help for mouse symbols...
### REF: https://github.com/waldronlab/HGNChelper/blob/master/vignettes/index.Rmd
### HGNChelper::checkGeneSymbols(stringr::str_to_title(df.genes_mouse_obesity$gene), unmapped.as.na = TRUE, map = NULL, species = "mouse") # Warning the list of valid mouse symbols seems to be incomplete, see below
df.genes_mouse_obesity
### Ortholog map
df.genes_mouse_obesity.human <- mouse_mgi_to_human_ortholog_mapping(df.genes_mouse_obesity, colname_genes_mgi="gene") 
df.genes_mouse_obesity.human <- df.genes_mouse_obesity.human %>% rename(ensembl_gene_id=ensembl_gene_id_human) # renaming column
head(df.genes_mouse_obesity.human)
### Write out
df.genes_mouse_obesity.human %>% write_tsv("yazdi2015_table1_mouse_obesity_genes.mapped.txt")

# 
# [1] "Reading gene mapping files..."
# [1] "Mapping from MGI to EnsemblID..."
# [1] "Converting MGI symbols to lowercase to avoid issues with case sensitivity"
# [1] "Mapping to orthologs..."
# [1] "Number of genes not mapped from MGI to Ensembl ID: 60 out of 224 genes (26.79 pct)"
# [1] "Number of genes not mapped to human ortholog: 65 out of 224 genes (29.02 pct)"
# [1] "Output dimension of final data table (ortholog mapped): n_genes=224 x n_features=6"
# [1] "List of genes not mapped from MGI to Ensembl ID:"
# [1] "adar2"     "alp1"      "asip"      "at2r"      "atx"       "cart"      "cb2r"      "chemr23"   "chop"      "cpt1"      "ctrp9"     "d2"        "dup(17)"  
# [14] "fatp4"     "fkbp51"    "foxo3a"    "gal3"      "girk4"     "gpr10"     "gpr120"    "gpr7"      "hif1α"     "hsd111b"   "hsd11β2"   "il18r"     "il1ri"    
# [27] "lrh1"      "magp1"     "mas"       "mt1a"      "mt1b"      "mthgh"     "nep"       "ngn3"      "p62"       "pc1"       "pc3"       "pcskin"    "pgds"     
# [40] "pi3k"      "pparg2"    "pparδ"     "ppir3a"    "ppkaa2"    "pref1"     "prrp"      "rage"      "ren"       "selm"      "sh2b"      "shp"       "spondin 2"
# [53] "src1"      "tap63"     "tis7"      "trkb"      "tw"        "xor"       "znt7"      "znt8"     
# [1] "List of genes not mapped to human ortholog:"
# [1] NA                   NA                   NA                   NA                   NA                   "ENSMUSG00000006464" NA                  
# [8] NA                   NA                   NA                   NA                   NA                   NA                   NA                  
# [15] NA                   NA                   NA                   NA                   "ENSMUSG00000020713" NA                   NA                  
# [22] NA                   NA                   NA                   NA                   NA                   NA                   NA                  
# [29] "ENSMUSG00000100916" NA                   NA                   NA                   NA                   NA                   NA                  
# [36] NA                   NA                   NA                   NA                   NA                   NA                   NA                  
# [43] NA                   "ENSMUSG00000032369" NA                   NA                   NA                   NA                   NA                  
# [50] "ENSMUSG00000021342" NA                   NA                   NA                   NA                   NA                   NA                  
# [57] NA                   NA                   NA                   NA                   NA                   NA                   NA                  
# [64] NA                   NA  




# ======================= file.genes_rare_variant ======================= #

### Human genes
file.genes_rare_variant <- "/projects/timshel/sc-genetics/sc-genetics/data/genes_obesity/original_tables/turcot2018_table1_rare_variants.txt"
df.genes_rare_variant <- read_tsv(file.genes_rare_variant)
df.genes_rare_variant.ens <- add_ensembl_gene_id_from_gene_symbol(df.genes_rare_variant, colname_genes="gene_symbol") # ---> all genes mapped
df.genes_rare_variant.ens %>% write_tsv("turcot2018_table1_rare_variants.mapped.txt")

# ======================= file.genes_mendelian ======================= #

file.genes_mendelian <- "/projects/timshel/sc-genetics/sc-genetics/data/genes_obesity/original_tables/turcot2018_s21_mendelian_obesity_genes.txt"
df.genes_mendelian <- read_tsv(file.genes_mendelian)
df.genes_mendelian
df.genes_mendelian.ens <- add_ensembl_gene_id_from_gene_symbol(df.genes_mendelian, colname_genes="gene_symbol") # ---> all genes mapped
df.genes_mendelian.ens %>% write_tsv("turcot2018_s21_mendelian_obesity_genes.mapped.txt")
