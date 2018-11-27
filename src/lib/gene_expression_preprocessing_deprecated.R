
library(Seurat)
library(Matrix)
library(parallel)


######################################################################################################
############################# Z-score related ################################
######################################################################################################

### FUNCTION WORKS AND TESTED.
wrapper_zscore_sems <- function(df.avg_expr) {
  ### INPUT: SEE calculate_and_write_sems
  print("Calculating Z-scores...")
  
  # ======= Z score: gene-wise z-score (each gene mean=0, sd=1) ======
  ### to give cell-type specificity
  ### to remove the effect of ubiquitous expressed genes.
  
  df.zscores.gene_wise <- as.data.frame(t(scale(t(as.data.frame(df.avg_expr)), center = TRUE, scale = TRUE))) # returns a DATA FRAME with columns=annotations, rownames=genes.
  df.zscores.gene_wise.na <- df.zscores.gene_wise %>% filter_all(any_vars(is.na(.))) # we (sometimes) have NA values. Don't know why. It could be because some genes have zero-variance
  if (nrow(df.zscores.gene_wise.na) > 0) {
    print(sprintf("OBS: N=%s genes removed because they contained NA values after zscore calculation. ", nrow(df.zscores.gene_wise.na)))
  }
  df.zscores.gene_wise <- df.zscores.gene_wise[complete.cases(df.zscores.gene_wise), ] # *OBS* remove genes with NA values
  
  ### pos_only
  df.zscores.gene_wise.pos_only <- replace(df.zscores.gene_wise, df.zscores.gene_wise < 0, 0) # replace negative values with 0
  ### log10(pos_only+1)
  df.zscores.gene_wise.pos_only_log <- log10(df.zscores.gene_wise.pos_only+1)
  write_csv(df.zscores.gene_wise.pos_only_log %>% rownames_to_column(var="gene"), sprintf("%s.z_score_pos_log.sem.csv.gz", name.dataset))
  
  # ======= Z score: cell-type-wise z-score (each cell-type mean=0, sd=1) ======
  ### to get the genes highest expressed within a cell-type.
  df.zscores.gene_wise.celltype_wise <- as.data.frame(scale(df.zscores.gene_wise)) # returns a DATA FRAME with columns=annotations, rownames=genes.
  write_csv(df.zscores.gene_wise.celltype_wise %>% rownames_to_column(var="gene"), sprintf("%s.z_score_two_step.sem.csv.gz", name.dataset))
  
  return(df.zscores.gene_wise)
}



### FULL CODE ZSCORE STEP-BY-STEP - WORKS. KEEP FOR DESCRIPTION OF RETURN VALUES
# df.avg_expr.t <- t(as.data.frame(df.avg_expr)) # transposing. returns a DATA FRAME with columns=genes, rownames=annotations
# df.avg_expr.t.scale <- scale(df.avg_expr.t, center = TRUE, scale = TRUE) # returns a *MATRIX* with colnumns=genes and rownames=annotations
# df.avg_expr.scale <- as.data.frame(t(df.avg_expr.t.scale)) # transposing back. returns a DATA FRAME with columns=annotations, rownames=genes.
# df.avg_expr.scale.na <- df.avg_expr.scale %>% filter_all(any_vars(is.na(.))) # we (sometimes) have NA values. Don't know why. It could be because some genes have zero-variance
# if (nrow(df.avg_expr.scale.na) > 0) {print(
#   sprintf("OBS: N=%s genes removed because they contained NA values after zscore calculation. ", nrow(df.avg_expr.scale.na)))
# }
# df.avg_expr.scale.clean <- df.avg_expr.scale[complete.cases(df.avg_expr.scale), ] # *OBS* remove genes with NA values
# write_csv(df.avg_expr.scale.clean %>% rownames_to_column(var="gene"), sprintf("%s.celltype_expr.z_score.hsapiens_orthologs.csv.gz", name.dataset))



######################################################################################################
############################################# SEM WRAPPERS ##############################################
######################################################################################################


# wrapper_mean_sems <- function(df.avg_expr, name.dataset, ortholog_mapping) {
#   ### INPUT
#   # df.avg_expr           data frame. genes x annotations. with rownamesand columnnames.
#   #                   as a 'standard' to ensure reproducibility, it is RECOMMENDED to *input a data.frame that has already been mapped to human genes*. 
#   # The various SEMs will change a bit depending on whether you do the mapping before or after calculating the SEMs. No one thing is more correct than the other, but best to do the same each time.
#   # name.dataset      prefix for output files
#   # ortholog_mapping  do automatic ortholog mapping.
#   #                   NULL = no mapping. Else a string to parse as 'type_mouse_gene_ids' argument to mouse_to_human_ortholog_gene_expression_mapping().
#   #                   valid strings "mgi" or "ensembl".
#   ### OUTPUT
#   # csv files written to working directory.
#   
#   print(sprintf("Ortholog mapping: %s", ortholog_mapping))
#   
#   # ======= Write Avg ======
#   write_sem(df.avg_expr, name.dataset, name.sem="avg_expr", ortholog_mapping)
#   
#   # ======= PEM ======
#   # df.pem <- calculate_pem(df.avg_expr)
#   # write_csv(df.pem %>% rownames_to_column(var="gene"), sprintf("%s.pem.sem.csv.gz", name.dataset))
#   
#   # ======= SI ======
#   df.SI <- specificity.index(df.avg_expr)
#   write_sem(df.avg_expr, name.dataset, name.sem="SI", ortholog_mapping)
#   
#   # ======= Z-scores ======
#   df.zscores.gene_wise <- wrapper_zscore_sems(df.avg_expr)
#   write_sem(df.avg_expr, name.dataset, name.sem="z_score", ortholog_mapping)
#   
#   # ======= Skene ======
#   df.specificity <- calculate_specificity(df.avg_expr)
#   write_sem(df.specificity, name.dataset, name.sem="specificity", ortholog_mapping)
#   
#   return(NULL)
# }


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


anno_specific_expr.avg_expr <- function(seurat_obj, colname_ident, genes=NULL, use.raw=F, do.log=T) {
  # DESCRIPTION: calculate average gene expression across cell types. 
  ### ***Seurat does the averaging in NON-log space*** ---> this is different from NEW November 2018 workflow. 
  
  # genes: gene names to calculate average on. If NULL, then Seurat will use all genes
  # If do.log_transform=T, ***log-transforms averaged data***. This CAN makes most sence, to dampen the gene expression variation. [NB the averaging in NON-log space by Seurat]
  # If use.raw=T, then the @raw.data slot is used instead of the CPM normalized @data slot
  # RETURN: data frame with columnnames=annotations, rownames=genes 
  
  print(sprintf("Running average expression for colname_ident=%s", colname_ident))
  
  ### AverageExpression(): by default uses `@data` slot (scaled=CPM_scaled and log-transformed single cell expression - NOT regressed), 
  # 'Output is in log-space when return.seurat = TRUE, otherwise it's in non-log space. Averaging is done in non-log space.'
  ### ***Seurat does the averaging in NON-log space***
  df.avg_expr <- AverageExpression(SetAllIdent(seurat_obj, id=colname_ident), genes.use=genes, use.raw=use.raw, return.seurat=F)
  # ---> RETURNS a DATA FRAME with columnnames=annotations, rownames=genes 
  
  if (do.log == TRUE) {
    print("Log transforming averaged data")
    df.avg_expr <- log1p(df.avg_expr)
  }
  return(df.avg_expr)
}


