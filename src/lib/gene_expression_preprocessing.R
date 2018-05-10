############### SYNOPSIS ###################
# Helper funtions for pre-processing expression data
# 1) annotation specific expression (e.g. cell type) 
# 2) ortholog gene mapping

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################

library(Seurat)
library(tidyverse)
library(Matrix)

library(parallel)


######################################################################################################
#################################### PRECOMPUTE T-STAT PARAMS ########################################
######################################################################################################
# --> write file
# cols: annotation, key=[other, annotation], value_var, value_mean, value_n



######################################################################################################
############################# Cell type specific expression functions ################################
######################################################################################################

TODO.FILTER_LOW_EXPRESSED_GENES <- function(TODO_NOT_IMPLEMENTED){
  print("TODO: write function to filter out low expressed genes prior to running t-test or avg. expression.")
}

anno_specific_expr.tstat.all_idents_parallel <- function(colname_ident, genes, seurat_obj, n.cores) {
  ### DESCRIPTION: runs anno_specific_expr.tstat.per_ident in parallel *across ALL idents in colname_ident*
  ### SEE anno_specific_expr.tstat.per_ident for details of arguments.
  ### OUTPUT
  # df.celltype_expr_ttest        data frame (genes x annotations) with t-statistics for annotation specific expression.
  #                               genes as rownames and annotations/cell-types as column name.
  
  
  start.time <- Sys.time() # time
  idents <- sort(unique(seurat_obj@meta.data[, colname_ident])) # get all idents
  print(sprintf("Running t-test for annotation specific expression for %s annotations", length(idents)))
  print(idents)
  print(sprintf("Will start cluster with n.cores=%s", n.cores))
  
  cl <- makeCluster(n.cores, type = "FORK")
  clusterEvalQ(cl, library(Seurat))
  clusterEvalQ(cl, library(tidyverse))
  list_of_dfs.celltype_expr_ttest <- parLapply(cl, idents, function(x_ident){anno_specific_expr.tstat.per_ident(colname_ident, x_ident, genes, seurat_obj)})
  stopCluster(cl)
  df.celltype_expr_ttest <- bind_cols(list_of_dfs.celltype_expr_ttest)
  # *OBS*: bind_cols(). combine list of dfs into a single data frame
  # ^ bind_cols(): rows are matched *by position*, so all data frames must have the same number of rows.
  rownames(df.celltype_expr_ttest) <- genes
  
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(sprintf("RUNTIME anno_specific_expr.tstat.all_idents_parallel: %.2f min (%.2f h)", time.taken/60, time.taken/(60*60) ))
  
  return(df.celltype_expr_ttest)
}

## # T-statistic cell type specific expression
anno_specific_expr.tstat.per_ident <- function(colname_ident, ident, genes, seurat_obj) {
  ### INPUT
  # colname_ident:                        name of any column in @meta.data that defines the grouping of ident (e.g. 'res.0.8' or 'annotation')
  # ident                                 a string or vector specifying the grouping. E.g. ident="beta" or ident=1
  # genes                                 character vector of gene names to perform t-test.
  ### Output
  # data frame                            single column data frame with t-statistics.
  #                                       the column name is the parameter value given by <ident> argument. Genes names are rowsnames.
  ### References
  # See implementation of Seurat DiffTTest : https://github.com/satijalab/seurat/blob/23853782a0793540a1db65cd12e95e90bfa86289/R/differential_expression.R#L1464
  # df.x <- DiffTTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1), cells.2 = WhichCells(object = pbmc_small, ident = 2)) # --> returns only p-values.
  
  ### For debugging
  #colname_ident <- "annotation"
  #ident <- "Astrocyte_1"
  
  print(sprintf("Running t-test cell-type specific expression for %s", ident))
  
  ### Get data
  data.test <- GetAssayData(object=seurat_obj, assay.type="RNA", slot="data")  # (same as seurat_obj@data?)
  # ^ returns a 'dgCMatrix'.
  # ^ rownames=gene names
  # ^ READ THIS: you can index using gene and cell names e.g. data.test["Gm1992", c("AAACCTGCAATCCGAT-1", "AAACCTGCAGATCCAT-1")], which returns a numeric vector that can be used for t-test
  
  ### Get cell labels
  seurat_obj <- SetAllIdent(seurat_obj, id=colname_ident) # set the ident
  cells.1 <- WhichCells(seurat_obj, ident=ident) # get the cells from ident [returns character vector of cell IDs]
  cells.2 <- WhichCells(seurat_obj, cells.use=setdiff(seurat_obj@cell.names, cells.1)) # the "rest" of the cells [returns character vector of cell IDs]
  
  ### Get all genes
  # genes <- rownames(data.test)[1:5] # TEST-CASE
  # genes <- rownames(data.test)
  
  ### Apply Welch two-sided t-test over genes
  # t_stats <- sapply(genes, function(gene) {t.test(x = data.test[gene, cells.1], y = data.test[gene, cells.2], alternative="two.sided", var.equal=F)$statistic}) # var.equal=F was used before March 29th 2018
  t_stats <- sapply(genes, function(gene) {t.test(x = data.test[gene, cells.1], y = data.test[gene, cells.2], alternative="two.sided", var.equal=T)$statistic}) 
  # lapply: lapply(X, FUN, ...) - returns a list (unnamed)
  # sapply: sapply(X, FUN, ...) - returns a vector (named with gene names)
  # TODO: consider using pblapply/pbsapply for status bar
  # (Welch --> var.equal=F)
  ### Create data frame and return
  df.res <- data.frame(t_stats, row.names=genes) # add rownames
  colnames(df.res) <- ident # add column name (e.g. "Astrocyte_1")
  return(df.res)
}


anno_specific_expr.avg_expr <- function(colname_ident, genes, seurat_obj) {
  # DESCRIPTION: calculate average gene expression across cell types. 
  # Returns average in log.space (this makes most sence, to dampen the gene expression variation)
  # RETURN: data frame with columnnames=annotations, rownames=genes 
  
  print(sprintf("Running average expression for colname_ident=%s", colname_ident))
  
  # AverageExpression(): uses `@data` slot (normalized and log-transformed single cell expression - NOT regressed)
  df.avg_expr <- log1p(AverageExpression(SetAllIdent(seurat_obj, id=colname_ident), genes.use=genes))
  # ---> RETURNS a DATA FRAME with columnnames=annotations, rownames=genes 
  
  return(df.avg_expr)
}










######################################################################################################
######################## mouse_to_human_ortholog_gene_expression_mapping #############################
######################################################################################################
# Funtion to do ortholog mapping: map mouse genes to human

### DESCRIPTION
# INPUT: avg expression file + mouse MGI-Ensembl mapping + ortholog mapping file.
# 1. Load avg expression file.
# 4. Map MOUSE GENE NAME to ENSEMBL: Map expression matrices MGI gene names to mouse Ensembl IDs
#     1. use timshel-lib from gene name to ensembl: "Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
# 5. SUBSET genes on 1-1 orthologs: in the expression file, keep only Ensembl IDs listed in "1-1 ortholog file"
# 6. REPLACE mouse ID with human ID.



mouse_to_human_ortholog_gene_expression_mapping <- function(df.expr) {
  ### INPUT
  # df.expr                 data frame containing expression data.
  #                         columns: cell-types or annotations
  #                         rows: genes (MGI symbols). The data frame *MUST have rownames* with *MGI symbols*.
  # file.out.map_stats:     output file path for the mapping stats (by default it just prints to the screen)
  ### OUTPUT
  # df.expr.ens.human:      a data frame where mouse genes have been substituted to human genes.
  #                         columns: cell-types or annotations
  #                         rows: genes (human Ensembl IDs) as rownames.

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
  # =========  EXPRESSION DATA: NUMERICAL TRANSFORMATION  ================= #
  # ======================================================================= #
  
  # NONE. We do this in the RolyPoly main script.
  
  # ======================================================================= #
  # ==========================  Map: MGI to EnsemblID  ==================== #
  # ======================================================================= #
  print("Mapping from MGI to EnsemblID...")
  
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
  str1 <- sprintf("Number of genes not mapped from MGI to Ensembl ID: %s out of %s genes (%.2f pct)", 
                  not_mapped_from_MGI_to_ensembl.n_genes, 
                  length(rownames(df.expr)),
                  not_mapped_from_MGI_to_ensembl.n_genes/length(rownames(df.expr))*100)
  str2 <- sprintf("Number of genes not mapped to human ortholog: %s out of %s genes (%.2f pct)", 
                  not_mapped_to_human_ortholog.n_genes,
                  length(rownames(df.expr.ens)),
                  not_mapped_to_human_ortholog.n_genes/length(rownames(df.expr.ens))*100)
  str3 <- sprintf("Output dimension of final data table (ortholog mapped): n_genes=%s x n_features=%s", 
                  nrow(df.expr.ens.human),
                  ncol(df.expr.ens.human))
  print(str1)
  print(str2)
  print(str3)
  print("List of genes not mapped from MGI to Ensembl ID:")
  print(not_mapped_from_MGI_to_ensembl.genes)
  print("List of genes not mapped to human ortholog:")
  print(not_mapped_to_human_ortholog.genes)
  
  return(df.expr.ens.human)
  
}  
  
  
  
