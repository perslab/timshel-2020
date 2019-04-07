############### SYNOPSIS ###################
# Helper funtions for gene mapping (orthologs, ensembl, entrz, ...)

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################


library(tidyverse)


######################################################################################################
######################## mouse_to_human_ortholog_gene_mapping [VECTOR INPUT] #############################
######################################################################################################

mouse_to_human_ortholog_gene_mapping <- function(gene_ids, type_mouse_gene_ids) {
  ### INPUT
  # gene_ids                ensembl or mgi mouse genes
  # type_mouse_gene_ids     "mgi" or "ensembl"
  ### OUTPUT
  # gene_ids.human:         a vector of same length as gene_ids but where mouse gene ids have been substituted to human gene ids.
  #                         mouse genes that did not map to human are assigned NA value.
  #                         hence you have do your mapping by ...
  if (!any(type_mouse_gene_ids %in% c("mgi", "ensembl"))) {
    stop(sprintf("Got wrong argument for type_mouse_gene_ids: %s", type_mouse_gene_ids))
  }
  print("Reading gene mapping files...")
  file.gene_name_mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
  df.gene_name_mapping <- suppressMessages(read_delim(file.gene_name_mapping, delim="\t"))
  file.ortholog_mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
  df.ortholog_mapping <- suppressMessages(read_delim(file.ortholog_mapping, delim="\t"))
  # ==========================  Map: MGI to EnsemblID  ==================== #
  # TODO: consider moving this step to a seperate function
  if (type_mouse_gene_ids=="mgi") {
    print("Mapping from MGI to EnsemblID...")
    print("Converting MGI symbols to lowercase to avoid issues with case sensitivity")
    df.gene_name_mapping$gene_name_optimal <- tolower(df.gene_name_mapping$gene_name_optimal)
    gene_ids <- tolower(gene_ids)
    ### Genes not found
    bool.not_found <- !(gene_ids %in% df.gene_name_mapping$gene_name_optimal) 
    not_mapped_from_MGI_to_ensembl.n_genes <- sum(bool.not_found) # n genes could not be mapped
    not_mapped_from_MGI_to_ensembl.genes <- gene_ids[bool.not_found] # e.g Fam150a not found, but this is because it is called ALKAL1. Fam150a is listed as a synonym in Emsembl now.
    ### Match
    idx.match <- match(gene_ids, df.gene_name_mapping$gene_name_optimal) # find mapping idx. 'nomatch' is set to NA
    gene_ids.ens <- df.gene_name_mapping$ensembl_gene_id[idx.match] # vector with same length as "gene_ids", but we NA values for genes not found.
  } else if (type_mouse_gene_ids=="ensembl") {
    gene_ids.ens <- gene_ids
  } else {
    stop("...")
  }
  # ========================  Subset on orthologs  ======================== #
  print("Subsetting on orthologs...")
  ### Find mouse genes with human orthologs
  idx.match <- match(gene_ids.ens, df.ortholog_mapping$mmusculus_homolog_ensembl_gene) # this works even when gene_ids.ens or $mmusculus_homolog_ensembl_gene have NA values (NA values are just 'probagated')
  gene_ids.ens.human <- df.ortholog_mapping$ensembl_gene_id[idx.match] # NA values for genes with no ortholog.
  ### Genes without orthologs
  bool.not_found <- !(gene_ids.ens %in% df.ortholog_mapping$mmusculus_homolog_ensembl_gene) 
  not_mapped_to_human_ortholog.n_genes <- sum(bool.not_found) # genes with no ortholog
  not_mapped_to_human_ortholog.genes <- gene_ids.ens[bool.not_found] # genes with no ortholog
  # ===============================  FINISH  ============================== #
  if (type_mouse_gene_ids=="mgi") {
    print("List of genes not mapped from MGI to Ensembl ID:")
    print(not_mapped_from_MGI_to_ensembl.genes)
  }
  print("List of genes not mapped to human ortholog:")
  print(not_mapped_to_human_ortholog.genes)
  if (type_mouse_gene_ids=="mgi") {
    str1 <- sprintf("Number of genes NOT mapped from MGI to Ensembl ID: %s out of %s genes (%.2f pct)", 
                    not_mapped_from_MGI_to_ensembl.n_genes, 
                    sum(!is.na(gene_ids)),
                    not_mapped_from_MGI_to_ensembl.n_genes/sum(!is.na(gene_ids))*100)
    print(str1)
  }
  str2 <- sprintf("Number of genes NOT mapped to human ortholog: %s out of %s genes (%.2f pct)", 
                  not_mapped_to_human_ortholog.n_genes,
                  sum(!is.na(gene_ids.ens)),
                  not_mapped_to_human_ortholog.n_genes/sum(!is.na(gene_ids.ens))*100)
  print(str2)
  str3 <- sprintf("Number of mapped genes in output: %s out of %s input genes (total mapping rate = %.2f pct)", 
                  sum(!is.na(gene_ids.ens.human)),
                  length(gene_ids.ens.human),
                  sum(!is.na(gene_ids.ens.human))/length(gene_ids.ens.human)*100)
  print(str3)
  return(gene_ids.ens.human)
}  


######################################################################################################
############ mouse_to_human_ortholog_gene_expression_mapping [DATA FRAME INPUT] ######################
######################################################################################################
# Funtion to do ortholog mapping: map mouse genes to human

### DESCRIPTION
# INPUT: avg expression file + mouse MGI-Ensembl mapping + ortholog mapping file.
# 1. Load avg expression file.
# 4. Map MOUSE GENE NAME to ENSEMBL: Map expression matrices MGI gene names to mouse Ensembl IDs
#     1. use timshel-lib from gene name to ensembl: "Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
# 5. SUBSET genes on 1-1 orthologs: in the expression file, keep only Ensembl IDs listed in "1-1 ortholog file"
# 6. REPLACE mouse ID with human ID.



mouse_to_human_ortholog_gene_expression_mapping <- function(df.expr, type_mouse_gene_ids="mgi") {
  ### INPUT
  # df.expr                 data frame containing expression data.
  #                         columns: cell-types or annotations
  # type_mouse_gene_ids     "mgi" or "ensembl"
  #                         rows: genes (MGI or mouse Ensembl symbols). The data frame *MUST have rownames* with *mouse genes*.
  # file.out.map_stats:     output file path for the mapping stats (by default it just prints to the screen)
  ### OUTPUT
  # df.expr.ens.human:      a data frame where mouse genes have been substituted to human genes.
  #                         columns: cell-types or annotations
  #                         rows: genes (human Ensembl IDs) as rownames.
  
  ### TODO
  # this function is VERY MEMORY INEFFICIENT because it copies many data frame. It is simple to fix.
  
  # ======================================================================= #
  # ==============================  INIT  ============================== #
  # ======================================================================= #
  
  if (!any(type_mouse_gene_ids %in% c("mgi", "ensembl"))) {
    stop(sprintf("Got wrong argument for type_mouse_gene_ids: %s", type_mouse_gene_ids))
  }
  
  ### make sure input is data frame (the below operations will fail on a tibble)
  df.expr <- as.data.frame(df.expr)
  
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
  if (type_mouse_gene_ids=="mgi") {
    print("Mapping from MGI to EnsemblID...")
    
    ### Converting MGI symbols to lowercase to avoid issues with case sensitivity
    print("Converting MGI symbols to lowercase to avoid issues with case sensitivity")
    df.gene_name_mapping$gene_name_optimal <- tolower(df.gene_name_mapping$gene_name_optimal)
    rownames(df.expr) <- tolower(rownames(df.expr))
    
    ### Genes not found
    bool.not_found <- !(rownames(df.expr) %in% df.gene_name_mapping$gene_name_optimal) 
    not_mapped_from_MGI_to_ensembl.n_genes <- sum(bool.not_found) # n genes could not be mapped
    not_mapped_from_MGI_to_ensembl.genes <- rownames(df.expr)[bool.not_found] # e.g Fam150a not found, but this is because it is called ALKAL1. Fam150a is listed as a synonym in Emsembl now.
    
    ### Match
    idx.match <- match(rownames(df.expr), df.gene_name_mapping$gene_name_optimal) # find mapping idx. 'nomatch' is set to NA
    genes.ensembl <- df.gene_name_mapping$ensembl_gene_id[idx.match] # vector with same length as "rownames(df.expr)", but we NA values for genes not found.
    
    
    ### Filter out genes not mapped, and replace gene names with IDs.
    df.expr.ens <- df.expr[!is.na(genes.ensembl),] # subset to mapped genes
    rownames(df.expr.ens) <- genes.ensembl[!is.na(genes.ensembl)] # replace gene names with IDs
  } else if (type_mouse_gene_ids=="ensembl") {
    df.expr.ens <- df.expr
  } else {
    stop()
  }
  
  # ======================================================================= #
  # ========================  Subset on orthologs  ======================== #
  # ======================================================================= #
  print("Subsetting on orthologs...")
  
  ### Find mouse genes with human orthologs
  idx.match <- match(rownames(df.expr.ens), df.ortholog_mapping$mmusculus_homolog_ensembl_gene)
  genes.human.orthologs <- df.ortholog_mapping$ensembl_gene_id[idx.match] # NA values for genes with no ortholog.
  
  ### Genes without orthologs
  bool.not_found <- !(rownames(df.expr.ens) %in% df.ortholog_mapping$mmusculus_homolog_ensembl_gene) 
  not_mapped_to_human_ortholog.n_genes <- sum(bool.not_found) # genes with no ortholog
  not_mapped_to_human_ortholog.genes <- rownames(df.expr.ens)[bool.not_found] # genes with no ortholog
  
  ### Subset genes
  df.expr.ens.human <- df.expr.ens[!is.na(genes.human.orthologs),] # subset to ortholog genes
  rownames(df.expr.ens.human) <- genes.human.orthologs[!is.na(genes.human.orthologs)] # replace mouse ensembl ID with human
  
  # ======================================================================= #
  # ===============================  FINISH  ============================== #
  # ======================================================================= #
  if (type_mouse_gene_ids=="mgi") {
    str1 <- sprintf("Number of genes not mapped from MGI to Ensembl ID: %s out of %s genes (%.2f pct)", 
                    not_mapped_from_MGI_to_ensembl.n_genes, 
                    length(rownames(df.expr)),
                    not_mapped_from_MGI_to_ensembl.n_genes/length(rownames(df.expr))*100)
    print(str1)
    print("List of genes not mapped from MGI to Ensembl ID:")
    print(not_mapped_from_MGI_to_ensembl.genes)
  }
  str2 <- sprintf("Number of genes not mapped to human ortholog: %s out of %s genes (%.2f pct)", 
                  not_mapped_to_human_ortholog.n_genes,
                  length(rownames(df.expr.ens)),
                  not_mapped_to_human_ortholog.n_genes/length(rownames(df.expr.ens))*100)
  print(str2)
  print("List of genes not mapped to human ortholog:")
  print(not_mapped_to_human_ortholog.genes)
  
  str3 <- sprintf("Output dimension of final data table (ortholog mapped): n_genes=%s x n_features=%s", 
                  nrow(df.expr.ens.human),
                  ncol(df.expr.ens.human))
  print(str3)
  
  
  
  return(df.expr.ens.human)
  
}  



######################################################################################################
############ mouse_to_human_ortholog_gene_expression_mapping [DATA FRAME INPUT] ######################
######################################################################################################
### THIS FUNCTION IS NOT USED VERY MUCH. [But it is still helpful.]

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




######################################################################################################
######################################## GENE ID MAPPING ############################################
######################################################################################################

# Map:  ENTREZ --> ENSEMBL
add_ensembl_ids_from_entrez <- function(df, colname_geneids_from="entrez", colname_geneids_to="ensembl_gene_id") {
  ### INPUT df: a tibble/data.frame with the column 'colname_geneids_from' with entrez gene ids.
  ### OUTOUT df 
  # returns a tibble with ensembl gene ids added to the column 'colname_geneids_to'. 
  # Genes that did not map have NA values in the 'colname_geneids_to' column.
  # If there are duplicated gene IDs in 'colname_geneids_to', then all but the first of the duplicated elements will be marked as NA.
  
  df <- as.data.frame(df) # convert to df to ensure the the below operations work.
  
  file.mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_id_mapping.hsapiens.ensembl_entrez.txt.gz"
  df.mapping <- suppressMessages(read_tsv(file.mapping))
  
  genes_mapped <- df.mapping$ensembl_gene_id[match(df[,colname_geneids_from], df.mapping$entrezgene)]
  bool_dups <- duplicated(genes_mapped, incomparables=NA) # marks elements with smaller subscripts as duplicates
  # ^ incomparables=NA: 'excluding' NA when counting duplicated. NA values will not be compared. (That is, duplicated() returns FALSE for NA values)
  # ^ duplicated(c(1,1,2,NA,NA,NA)) returns FALSE  TRUE FALSE FALSE  TRUE  TRUE.
  # ^ duplicated(c(1,1,2,NA,NA,NA), incomparables=NA) returns FALSE  TRUE FALSE FALSE FALSE FALSE.
  print(sprintf("Number of genes mapped: %s",sum(!is.na(genes_mapped))))
  print(sprintf("Number of genes not mapped: %s",sum(is.na(genes_mapped)))) # number of not mapped genes
  print(sprintf("Number of genes with a NON-unique mapping (genes with duplicated ensembl gene IDs after mapping): %s",sum(bool_dups)))
  ### set duplicated rows (with smaller subscripts) as NA
  genes_mapped[bool_dups] <- NA
  print(sprintf("Total mapping stats: %s genes have no mapping (not mapped + duplicates) out of %s input genes.", sum(is.na(genes_mapped)), length(genes_mapped)))
  print(sprintf("Total genes mapped (non NA genes): %s", sum(!is.na(genes_mapped))))
  df <- df %>% mutate(!!rlang::sym(colname_geneids_to):=genes_mapped) %>% as.tibble()
  # filter(!is.na(gene)) %>% # remove all rows without mapping
  # filter(!duplicated(gene)) # keep only one of the duplicated pair (if any)
  print(sprintf("Returning tibble with the column '%s' added where all gene identifiers unique. Unmapped genes have NA values", colname_geneids_to))
  return(df)
}


# Map: HUMAN ENSEMBL --> GENE SYMBOL
hs_add_gene_symbol_from_ensembl_ids <- function(df, colname_geneids_from="ensembl_gene_id", colname_geneids_to="gene_symbol") {
  ### INPUT df: a tibble/data.frame with the column 'colname_geneids_from' with human ensembl gene ids.
  ### OUTOUT df 
  # returns a tibble with human gene symbols added to the column 'colname_geneids_to'. 
  # Genes that did not map have NA values in the 'colname_geneids_to' column.
  # If there are duplicated gene IDs in 'colname_geneids_to', then all but the first of the duplicated elements will be marked as NA.
  
  df <- as.data.frame(df) # convert to df to ensure the the below operations work.
  
  file.mapping <- "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/GRCh38.ens_v90.gene_name_version2ensembl.txt.gz"
  df.mapping <- suppressMessages(read_tsv(file.mapping))
  
  genes_mapped <- df.mapping$gene_name_optimal[match(df[,colname_geneids_from], df.mapping$ensembl_gene_id)]
  bool_dups <- duplicated(genes_mapped, incomparables=NA) # marks elements with smaller subscripts as duplicates
  # ^ incomparables=NA: 'excluding' NA when counting duplicated. NA values will not be compared. (That is, duplicated() returns FALSE for NA values)
  # ^ duplicated(c(1,1,2,NA,NA,NA)) returns FALSE  TRUE FALSE FALSE  TRUE  TRUE.
  # ^ duplicated(c(1,1,2,NA,NA,NA), incomparables=NA) returns FALSE  TRUE FALSE FALSE FALSE FALSE.
  print(sprintf("Number of genes mapped: %s",sum(!is.na(genes_mapped))))
  print(sprintf("Number of genes not mapped: %s",sum(is.na(genes_mapped)))) # number of not mapped genes
  print(sprintf("Number of genes with a NON-unique mapping (genes with duplicated ensembl gene IDs after mapping): %s",sum(bool_dups)))
  ### set duplicated rows (with smaller subscripts) as NA
  genes_mapped[bool_dups] <- NA
  print(sprintf("Total mapping stats: %s genes have no mapping (not mapped + duplicates) out of %s input genes.", sum(is.na(genes_mapped)), length(genes_mapped)))
  print(sprintf("Total genes mapped (non NA genes): %s", sum(!is.na(genes_mapped))))
  df <- df %>% mutate(!!rlang::sym(colname_geneids_to):=genes_mapped) %>% as.tibble()
  # filter(!is.na(gene)) %>% # remove all rows without mapping
  # filter(!duplicated(gene)) # keep only one of the duplicated pair (if any)
  print(sprintf("Returning tibble with the column '%s' added where all gene identifiers unique. Unmapped genes have NA values", colname_geneids_to))
  return(df)
}

