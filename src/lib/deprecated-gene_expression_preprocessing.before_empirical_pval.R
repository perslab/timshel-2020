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


library(tidyverse)


######################################################################################################
######################################## Binning ############################################
######################################################################################################

wrapper_binning <- function(df, n_bins=101, threshold_bin_zero=0){
  print(sprintf("Binning data frame with n_bins=%s", n_bins))
  df.binned <- df %>% transmute_all(funs(bin_vector), n_bins, threshold_bin_zero) # *OBS* output tibble looses rownames
  # rownames(df.binned) <- rownames(df) # set rownmaes *IMPORTANT*
  return(df.binned)
}

bin_vector <- function(vector_in, n_bins, threshold_bin_zero){
  # Function to compute bins 
  # INPUT: a vector of numeric values to bin
  # OUTPUT: a numeric vector with values {0 ... n_bins-1?}
  # Skene used n_bins=41
  # *REMARK #1*: note that vector_in <= threshold_bin_zero will be put in bin 0.
  # *REMARK #2*: note that NA values in vector_in will be put in bin 0.
  
  bin_values <- rep(0,length(vector_in))
  bool.not_bin_zero <- (vector_in > threshold_bin_zero) & (!is.na(vector_in)) # bool index for elements that should not go in bin 0
  bin_values[bool.not_bin_zero] <- as.numeric(cut(vector_in[bool.not_bin_zero], # cut divides the range of x into intervals and codes the values in x according to which interval they fall. 
                                                                 breaks=unique(quantile(vector_in[bool.not_bin_zero], probs=seq(0,1, by=1/n_bins))), 
                                                                 include.lowest=TRUE))
  ### ORIGINAL CODE (by Skene)
  # quantileValues = rep(0,length(vector_in))
  # quantileValues[vector_in>0] = as.numeric(cut(vector_in[vector_in>0], # cut divides the range of x into intervals and codes the values in x according to which interval they fall. 
  #                                              breaks=unique(quantile(vector_in[vector_in>0], probs=seq(0,1, by=1/numberOfBins), na.rm=TRUE)), 
  #                                              include.lowest=TRUE))
  
  return(bin_values)
}


######################################################################################################
############################# SKENE expression functions ################################
######################################################################################################


# wrapper_specificity_quantiles <- function(df.avg_expr, name.dataset) {
#   ### INPUT: data frame or matrix of average expression. genes x celltypes. with rownames and columnnames.
#   # SEE calculate_specificity() for notes on sensitivity to log_transformation.
#   ### OUTPUT: data frame with specificity quantiles.
#   
#   df.specificity <- calculate_specificity(df.avg_expr)
#   write_csv(df.specificity %>% rownames_to_column(var="gene"), sprintf("%s.specificity.sem.csv.gz", name.dataset))
#   
#   df.specificity.quantiles <- df.specificity %>% transmute_all(funs(transform_to_quantiles), threshold_bin_zero=0) # *OBS* output tibble looses rownames
#   rownames(df.specificity.quantiles) <- rownames(df.specificity) # set rownmaes *IMPORTANT*
#   write_csv(df.specificity.quantiles %>% rownames_to_column(var="gene"), sprintf("%s.specificity_quantiles.sem.csv.gz", name.dataset))
# }


calculate_specificity <- function(df.avg_expr){
  ### INPUT:
  # df.avg_expr:      data frame or matrix of average expression. genes x celltypes. with rownames and columnnames.
  #                   It is NOT recommended to input  average expression in log space, but the specificity measure is robust to log-transformation, so it will only make very small difference.
  ### OUTPUT
  # df.specificity     data.frame. genes x cell_types. with rownames and columnnames
  #                     each row (gene) in the matrix sums to 1.
  ### REF 
  # EWCE package: calculate.specificity.for.level() inside script generate.celltype.data.r
  # https://github.com/NathanSkene/EWCE/blob/47c429c2e5983d7a6a8ba70e7f644572d4c2ef72/R/generate.celltype.data.r
  ### EXPLANATION
  # Skenes method is "dividing the mean UMI counts for a gene A in cell type X with the mean UMI counts for gene A for all cells". 
  # Hence the specificity measure is between 0 and 1. 
  # The measure is expressing "what percentage of gene A's expression is coming from cell type X". 
  # So if all cell types express gene A equally, all the cell-types will have more or less the same specificity measure for that gene.
  print("Calculating specificity (Skene)")
  
  ### Scaling every cells expression, by the sum of the cells expression.
  ### This is basically like setting scale.factor=1 in Seurat NormalizeData()
  mat.mean_expr.scaled = t(df.avg_expr)/colSums(df.avg_expr) # celltypes x genes.  [Write out the matrices if needed for understanding]
  mat.mean_expr.scaled <- t(mat.mean_expr.scaled) # genes x celltypes. Every column should now sum to 1.
  # ^ I think the above is equivalent to sweep(df, 2, colSums(df), `/`)
  
  ### Calculate specificity
  df.specificity = as.data.frame(mat.mean_expr.scaled/(rowSums(mat.mean_expr.scaled)+0.000000000001))
  # ^ dividing every 'cell expression profile column' with the genes total expression (the rowSums)
  # ^ adding 0.000000000001 avoids 'devision by zero errors'
  # ^NB: rowSums(mat) == apply(mat,1,sum)
  return(df.specificity)
  
  
  ### ORIGINAL SKENE FUNCTION
  # calculate.specificity.for.level <- function(df.mean_exp){
  #   normalised_meanExp = t(t(df.mean_exp)*(1/colSums(df.mean_exp)))
  #   specificity = normalised_meanExp/(apply(normalised_meanExp,1,sum)+0.000000000001)
  #   return(specificity)
  # }
}




######################################################################################################
############################################# SI ORIGINAL ##############################################
######################################################################################################

normalized_specificity_index <- function (df) {
  ### Input: data.frame. genes x cell-types. Any rownames are not used
  ### Output: normalized SI data.frame with same dimensions as input df.
  ### Algorithm inspired by specificity.index() function from pSI packages: http://genetics.wustl.edu/jdlab/psi_package/
  ### Improvements of nSI (normalized specificity index) compared to the original SI
  # Summary: modified for single-cell data.
  # 1) Normalized ranks: if a gene is only expressed in one cell-type it will get nSI = 1
  # 2) Added epsilon to ratio: important for de
  
  print("Calculating SI")
  epsilon <- 1e-100
  df <- as.data.frame(df) # make sure the input is a data frame.
  n_annotations <- ncol(df) 
  df.SI <- data.frame(matrix(nrow = nrow(df) , ncol = n_annotations)) 
  colnames(df.SI) <- colnames(df) 
  # rownames(df.SI) <- rownames(df) # PT note: rownames are used
  for (j in 1:n_annotations) { # j: index cell-types (columns of input)
    print(sprintf("Calculating SI: #%s/#%s", j, n_annotations))
    notme <- c(1:n_annotations)[-j] # exclude self index
    ### Calculate fold-change
    # fc <- df[, j]/df[, notme] # without epsilon
    fc <- (df[, j]+epsilon)/(df[, notme]+epsilon) # df, genes x (annotations-1) | column-wise division.
    # ^ we add epsilon because if a gene is not expressed in the other annotations, the ratio will be NA/Inf.
    # ^ then when we rank the genes while ignoring NA values (na.last="keep"), the gene will get a low rank or maybe even be NA.
    # ^ fc=~0 --> 0: ok
    # ^ fc=~1 --> 0: ok, see below
    # ^ fc>>c --> do nothing: fc can get very big (e.g. ~1e100) when annotation_j > 0 and annotation_k is 0: (some_number+epsilon)/(0+epsilon). We don't need to do anything there. The ranks will still work.
    fc[dplyr::near(fc,0)] <- 0 # we map fc near 0 to 0 to avoid numerical precision artifacts when ranking the genes and keeping track of ties. fc values near 0 arise when the annotation[j] is zero and the epsilon is added: (0+epsilon)/(some_number+epsilon)
    fc[dplyr::near(fc,1)] <- 0 # we map fc near 1 to 0 because these fc values are almost certainly a result of dividing two zero measurements: (0+epsilon)/(0+epsilon). These fc values should be 0, because the gene is not expressed in the given annotation.
    ### Rank genes: would want high ranks to correspond to high specificity
    fc_ranks <- apply(fc, 2, base::rank, na.last="keep", ties.method="min") # matrix, (genes x annotations-1)
    # ^ rank in ascending order (low values = low rank; high values = high rank)
    # ^ na.last="keep": NA values are kept with rank NA. That is, NA values are ignored during ranking and 'propagated' to the result vector.
    # ^ ties.method="min" : tied genes get the minimum rank. We do this to add maximal penalty to genes with fc = 0. We are awere that it will also effect genes with fc >>, but they are unlikely to have tied values 
    # ALT FASTER: fc_ranks <- apply(fc, 2, data.table::frank, ...) # identical but faster solution using data.table::frank
    ### Normalize ranks: divide each column by the number of genes
    fc_ranks_normalized <- rowMeans(fc_ranks/nrow(fc_ranks))
    df.SI[,j] <- fc_ranks_normalized 
  }
  return(df.SI)
}


specificity.index.original <- function (pSI.in) {
  ### Input: data.frame. genes x cell-types. Must include row names and colnames.
  ### returns only SI data.frame. with same rownames and colnames as input df.
  ### Modified copy of specificity.index() function from pSI packages.
  print("Calculating SI")
  
  pSI.in <- as.data.frame(pSI.in) # make sure the input is a data frame.
  
  lng <- length(pSI.in[1, ])
  datComb <- array(NA, c(length(pSI.in[, 1]), 2 * lng))
  colnames(datComb) <- c(rep(colnames(pSI.in), each = 2))
  psi_columns <- c(seq(from = 2, to = ncol(datComb), by = 2))
  si_columns <- c(seq(from = 1, to = ncol(datComb), by = 2))
  rownames(datComb) <- rownames(pSI.in) # PT note: rownames are used
  for (j in 1:lng) { # j: index cell-types (columns of input)
    print(sprintf("Calculating SI: #%s/#%s", j, lng))
    notme = c(1:(j - 1), (j + 1):lng)
    if (j == 1) {
      notme = c(2:lng)
    }
    if (j == lng) {
      notme = c(1:(lng - 1))
    }
    
    zO <- array(NA, c(length(pSI.in[, 1]), lng - 1))
    zO <- log2(pSI.in[, j]/pSI.in[, notme])
    rankzO <- array(0, dim(zO))
    zO <- zO * -1 # to invert the data, so low ranks end up being 'good', and NAs will go to the bottom?
    for (i in 1:(lng - 1)) {
      rankzO[, i] = rank(zO[, i], na.last = "keep")
    }
    avg_ipO <- c()
    avg_ipO <- .rowMeans(rankzO, m = nrow(rankzO), n = ncol(rankzO))
    datComb[, ((j * 2) - 1)] <- avg_ipO # SI?
  } # END LOOP 'for (j in 1:lng)'. Index cell-types (columns of input)
  df.SI <- data.frame(datComb[, si_columns])
  return(df.SI)
}


######################################################################################################
############################################# PEM ##############################################
######################################################################################################

# Preferential Expression Measure
# * PEM estimates how different the expression of the gene is relative to an expected expression, under the assumption of uniform expression in all tissues. 
# * Formula: PEM(gene, c_i) = log10(C/c_i * g/G_i)
# * REF derivation: https://academic.oup.com/bib/article/18/2/205/2562739
# * REF orig: https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-4-31

# * NB all data is averaged across cell-types, so each cell-type only have one column.
# * C: sum expression of all genes in all cell-types. [CONSTANT]
# * c_i:  sum expression of all genes in cell-type i (self).
# * G_i: sum expression of gene in all cell-types.
# * g: expression of gene in cell-type i (self).

# * Pseudo code
# * df : genes x cell-types. Add rownames.
# * C = sum(sum(df)) # total sum
# * c = colSums(df) # celltype sum | 1 x cell-types. Named for indexing.
# * G = rowSums(df) # gene sum | 1 x genes. Named for indexing.
# * foreach gene j:
#   * foreach celltype i:
#   * PEM = log10(C/c[i]* df[i,j]/G[j]

calculate_pem <- function(df.avg_expr){
  print("Calculating PEM...")
  total_sum <- sum(sum(df.avg_expr)) # total sum
  celltype_sums <- colSums(df.avg_expr) # celltype sum | 1 x cell-types. Named for indexing.
  gene_sums <- rowSums(df.avg_expr) # gene sum | 1 x genes. Named for indexing.
  df.pem <- data.frame(matrix(ncol = ncol(df.avg_expr), nrow=nrow(df.avg_expr)), row.names = rownames(df.avg_expr))
  colnames(df.pem) <- colnames(df.avg_expr)
  counter <- 1
  for (celltype in colnames(df.pem)) {
    print(sprintf("Calculating PEM: #%s/#%s", counter, length(colnames(df.pem))))
    for (gene in rownames(df.pem)) {
      df.pem[gene, celltype] = log10(total_sum/celltype_sum[celltype]* df.avg_expr[gene,celltype]/gene_sums[gene]) # add +1 to avoid -inf for log(0)
    }
    counter <- counter + 1
  }
  return(df.pem)
}



######################################################################################################
################################## SEM CALC that takes SEM as input  ##################################
######################################################################################################

zscore_sem <- function(sem) {
  print("Calculating z-score")
  df.zscores.gene_wise <- as.data.frame(t(scale(t(as.data.frame(sem[["mean"]])), center = TRUE, scale = TRUE))) # returns a DATA FRAME with columns=annotations, rownames=genes.
  df.zscores.gene_wise.na <- df.zscores.gene_wise %>% filter_all(any_vars(is.na(.))) # we (sometimes) have NA values. Don't know why. It could be because some genes have zero-variance
  if (nrow(df.zscores.gene_wise.na) > 0) {
    print(sprintf("OBS: N=%s genes contained NA values after zscore calculation. ", nrow(df.zscores.gene_wise.na)))
  }
  return(df.zscores.gene_wise)
}

tstat_sem <- function(sem) {
  ### vectorized calculation of pooled variance for each gene
  # formula: sum_i{var[i]*(n[i]-1)}/sum_i{n[i]-1}, where i is the group.
  print("Calculating t-stat")
  df.n <- sem[["ncells"]] %>% slice(rep(1:n(), each=length(sem[["genes"]]) )) # df, 1 x annotations | replicate row in ncells. REF rep(): https://stackoverflow.com/a/47780537/6639640
  var.pooled <- rowSums(sem[["var"]] * (df.n-1))/rowSums(df.n - 1) # numeric vector, genes x 1 | unbiased least squares estimate (n-1) of pooled variance.
  list.tstat <- list()
  for (annotation in sem[["annotations"]]) {
    print(annotation)
    annotation_sym <- rlang::sym(annotation)
    x <- sem[["mean"]] %>% pull(!!annotation_sym) # numeric vector, genes x 1 | mean self
    ### Calculate mean for 'other' using ncells as weights to recover the true mean of the 'other' group.
    # formula: sum_i{mean[i]*n[i]}/sum_i{n[i]}
    df.n_other <- sem[["ncells"]] %>% select(-!!annotation_sym) %>% slice(rep(1:n(), each=length(sem[["genes"]]) )) # df, genes x annotations. Each row contains the number of cells and each row is identical.
    y <- rowSums(sem[["mean"]] %>% select(-!!annotation_sym) * df.n_other)/rowSums(df.n) # numeric vector, genes x 1
    # ^ df.mean_other * df.n_other: column-wise multiplication. 
    # ^ rowSums(df.n %>% slice(1)) gives the same as rowSums(df.n) because all rows are identical
    tstat <- (x-y)/var.pooled # numeric, genes x 1
    # pt(tstat, df=sum(df.n)-2, lower.tail=T) # numeric vector, genes x 1 | one-sided ttest for higher expression. Each gene as the same number of degree of freedom (dof). Dof is total number of observations (cells) minus 2.
    list.tstat[[annotation]] <- tstat
    # if (annotation=="ENT1") {break}
  }
  df.tstat <- bind_cols(list.tstat) # df/tbl, genes x annotation | binds by position
  return(df.tstat)
}



ges_sem <- function(sem) {
  ### Calculate Gene Enrichment Score
  ### Ref: Zeisel, Cell 2018
  ### formula: (mean[celltype]+epsilon_1)/(mean[other]+epsilon_1) * (frac[celltype]+epsilon_2)/(frac[other]+epsilon_2)
  ### epsilon_1=0.1 ; epsilon_2=0.01
  print("Calculating GES")
  epsilon_1 <-  1e-100 # or .Machine$double.eps
  epsilon_2 <- 0.01
  ### vectorized calculation of pooled variance for each gene
  # formula: sum_i{var[i]*(n[i]-1)}/sum_i{n[i]-1}, where i is the group.
  df.n <- sem[["ncells"]] %>% slice(rep(1:n(), each=length(sem[["genes"]]) )) # df, 1 x annotations | replicate row in ncells. REF rep(): https://stackoverflow.com/a/47780537/6639640
  list.ges <- list()
  for (annotation in sem[["annotations"]]) {
    print(annotation)
    annotation_sym <- rlang::sym(annotation)
    mean.self <- sem[["mean"]] %>% pull(!!annotation_sym) # numeric vector, genes x 1
    frac.self <- sem[["frac_expr"]] %>% pull(!!annotation_sym) # numeric vector, genes x 1
    ### Calculate mean and frac for 'other' using ncells as weights
    df.n_other <- sem[["ncells"]] %>% select(-!!annotation_sym) %>% slice(rep(1:n(), each=length(sem[["genes"]]) )) # df, genes x annotations. Each row contains the number of cells and each row is identical.
    mean.other <- rowSums(sem[["mean"]] %>% select(-!!annotation_sym) * df.n_other)/rowSums(df.n) # numeric vector, genes x 1
    frac.other <- rowSums(sem[["frac_expr"]] %>% select(-!!annotation_sym) * df.n_other)/rowSums(df.n) # numeric vector, genes x 1
    ### Calculate GES
    ges <- (mean.self+epsilon_1)/(mean.other+epsilon_1) * (frac.other+epsilon_2)/(frac.other+epsilon_2) # numeric, genes x 1
    ges[dplyr::near(ges, 0)] <- 0 # assigning ges values 'near equal to' zero to zero.
    list.ges[[annotation]] <- ges
  }
  df.ges <- bind_cols(list.ges) # df/tbl, genes x annotation | binds by position
  return(df.ges)
}


######################################################################################################
############################################# SEM UTILS ##############################################
######################################################################################################


read_file_fast <- function(file_path) {
  print(sprintf("Reading file %s", file_path))
  # fast reading of files
  df <- data.table::fread(file_path, nThread=24, showProgress=T) %>% as.tibble()
  # file_path can be a gzip file. It will automatically be infered: https://github.com/hansenlab/bsseq/commit/0656078974325effacdaabeaba68bc5c07bfe704
  # Compressed files ending .gz and .bz2 are supported if the R.utils package is installed.
  return(df)
}


write_sems <- function(object, slot, name.dataset) {
  ### Function: write SEM: "sem" or "sem_bin"
  ### We recommend that name.dataset contains information about specie
  stopifnot(slot %in% c("sem", "sem_bin"))
  print(sprintf("Writing slot=%s files for name.dataset=%s", slot, name.dataset))
  for (name.sem in names(object[[slot]])) {
    file.out <- sprintf("%s.%s.csv.gz", name.dataset, name.sem)
    print(sprintf("Writing file: %s", file.out))
    object[[slot]][[name.sem]] %>% mutate(gene=object[["genes"]]) %>% select(gene, everything()) %>% write_csv(path=file.out)
    # or use data.table::fwrite() and afterwards R.utils::gzip('filename.csv',destname='filename.csv.gz')
  }
}

### OLD
# write_sem <- function(df, name.dataset, name.sem, ortholog_mapping=NULL) {
#   print(sprintf("Writing files for SEM = %s", name.sem))
#   write_csv(df %>% rownames_to_column(var="gene"), sprintf("%s.%s.sem.csv.gz", name.dataset, name.sem))
#   
#   if (!is.null(ortholog_mapping)) {
#     df.human <- mouse_to_human_ortholog_gene_expression_mapping(df, type_mouse_gene_ids=ortholog_mapping)
#     write_csv(df.human %>% rownames_to_column(var="gene"), sprintf("%s.%s.hsapiens_orthologs.sem.csv.gz", name.dataset, name.sem))
#   }
# }

identical_value <- function(x,y) {
  # returns NA if any values are not identical
  if (identical(x,y)) x else NA # or x==y instead of identical(). 
}

######################################################################################################
############################################# SEM_OBJECT METHODS ##############################################
######################################################################################################

### sem_object
# var
# frac_expr
# mean
# genes_exclude
# ncells
# genes --> character
# annotations --> character
# sem --> list(sems)
# sem_bin --> list(sems)
# group_by_annotation.sem_bin --> list of dfs (annotations)
# group_by_annotation.sem --> list of dfs (annotations)
# sem_meta --> list of dfs (mean,median,sd)  made from group_by_annotation.sem_bin


preprocess_and_set_slots <- function(object) {
  ### check that genes and annotations are identical
  ### run only once because 'gene' column is removed.
  list_annotations <- list()
  list_genes <- list()
  for (x in c("frac_expr", "mean", "var")) {
    if (anyNA(object[[x]])) {
      stop(sprintf("Input data contains NA values: %s", x))
    }
    list_genes[[x]] <- object[[x]] %>% pull(gene)
    object[[x]] <- object[[x]] %>% select(-gene) # remove col
    list_annotations[[x]] <- colnames(object[[x]]) # set annotations after removing gene col
  }
  if (anyNA(Reduce(f=identical_value, x=list_annotations))) {
    stop("Annotations are not identical for pre-loaded data. Input data must contain identical cell-types in the same order.")
  }
  if (anyNA(Reduce(f=identical_value, x=list_genes))) {
    stop("Genes are not identical for pre-loaded data. Input data must contain identical genes in the same order.")
  }
  
  ### Set genes and annotations
  object[["genes"]] <- list_genes[[1]]
  object[["annotations"]] <- list_annotations[[1]]
  
  ### Last check
  stopifnot(all(object[["annotations"]] == colnames(object[["ncells"]]))) # annotation must match annotations in ncell
  return(object)
}



create_sem_object <- function(file_prefix) {
  object <- list()
  class(object) <- "sem object"
  object[["var"]] <- read_file_fast(sprintf("%s.pre_calc.var.csv.gz", file_prefix))
  object[["frac_expr"]] <- read_file_fast(sprintf("%s.pre_calc.frac_expr.csv.gz", file_prefix))
  object[["mean"]] <- read_file_fast(sprintf("%s.pre_calc.mean.csv.gz", file_prefix))
  object[["genes_exclude"]] <- read_file_fast(sprintf("%s.pre_calc.sporadically_expressed_genes.anova.csv.gz", file_prefix))
  object[["ncells"]] <- read_file_fast(sprintf("%s.pre_calc.ncells.csv.gz", file_prefix)) %>% 
    rename(annotation=cell_type, n=NCells) %>% # *TMP*
    distinct() # *TMP* # df, 2 x annotations | cols = {annotation, n}
  
  ### making 'ncells' a one row df.
  object[["ncells"]] <- object[["ncells"]] %>% 
    spread(key=annotation, value=n) %>% # 1-row tibble with annotation as colnames. spread() will re-arrange column order, so we need to fix that in the next step
    select(object[["ncells"]] %>% pull(annotation)) # df, annotations x 1 | re-arrrenge columns back to their original order.
  ### alternative WORKS.
  # n <- object[["ncells"]] %>% pull(n)
  # names(n) <- object[["ncells"]] %>% pull(annotation)
  # object[["ncells"]] <- as.tibble(as.list(n))
  
  ### create object slots (not needed)
  # object[["sem"]] <- list()
  # object[["sem_bin"]] <- list()
  
  ### pre_process
  object <- preprocess_and_set_slots(object) 
  
  return(object)
}


subset_annotations <- function(object, annotations) {
  ### Function: returns a object with the subset of annotations provided in 'annotations'.
  ### annotations: a character vector of annotation names to subset on. 
  ### NB: annotations are re-ordered by the order provided in annotations.
  print(sprintf("Subsetting object: keeping n=%s annotations.", length(annotations)))
  
  list_annotations <- list()
  for (x in c("frac_expr", "mean", "var", "ncells")) {
    object[[x]] <- object[[x]] %>% select(!!!rlang::syms(annotations)) # Notice !!!syms() because annotations is a vector.
    list_annotations[[x]] <- colnames(object[[x]])
  }
  if (anyNA(Reduce(f=identical_value, x=list_annotations))) {
    stop("Annotations are not identical for pre-loaded data. Input data must contain identical cell-types in the same order.")
  }
  object[["annotations"]] <- object[["annotations"]][object[["annotations"]] %in% annotations] # subject annoations
  
  ### Excluding zero mean genes (after subsetting on annotations)
  bool.genes_exclude <- rowMeans(object[["mean"]])==0
  print(sprintf("Number of zero-mean genes marked for exclusion: %s", sum(bool.genes_exclude)))
  for (x in c("frac_expr", "mean", "var")) {
    object[[x]] <- object[[x]] %>% filter(!bool.genes_exclude)
  }
  object[["genes"]] <- object[["genes"]][!bool.genes_exclude]
  print(sprintf("Number of genes in object: %s", length(object[["genes"]])))
  
  ### Clear data slots
  print("Clearning sem data slots: sem, sem_bin, etc...")
  object[["sem"]] <- NULL
  object[["sem_bin"]] <- NULL
  object[["group_by_annotation.sem"]] <- NULL
  object[["group_by_annotation.sem_bin"]] <- NULL
  object[["sem_meta"]] <- NULL
  return(object)
}

exclude_genes <- function(object) {
  ### exclude the following genes
  ### 1) genes with insignificant anova test
  ### 2) genes with mean zero expression in ALL annotation groups.
  ### OBS: some genes in 'genes_exclude' have 'NA' values in the pvalue column.
  ###       this is because ALL annotation groups have zero variance for that gene.
  ###       these genes will almost always be genes with zero mean expression.
  ###
  genes.anova_na <- object[["genes_exclude"]] %>% filter(is.na(pvalue)) %>% pull(gene)
  print(sprintf("Number of genes with NA in anova pvalue: %s. (This number is not used for anything. We just report it for your information)", length(genes.anova_na))) # 
  
  ### Anova genes
  threshold.sporadically_expressed <- 0.00001 # Skene threshold 0.00001
  genes.sporadic_expressed <- object[["genes_exclude"]] %>% filter(!is.na(pvalue) & (pvalue > threshold.sporadically_expressed)) %>% pull(gene)
  print(sprintf("Number of sporadic_expressed_genes: %s", length(genes.sporadic_expressed)))
  bool.sporadic_expressed <- object[["genes"]] %in% genes.sporadic_expressed
  
  ### Zero mean genes
  bool.genes_zero_mean <- rowMeans(object[["mean"]])==0
  print(sprintf("Number of zero-mean genes: %s", sum(bool.genes_zero_mean)))
  
  bool.genes_exclude <- bool.genes_zero_mean | bool.sporadic_expressed 
  print(sprintf("Total number of genes marked for exclusion: %s", sum(bool.genes_exclude)))
  for (x in c("frac_expr", "mean", "var")) {
    object[[x]] <- object[[x]] %>% filter(!bool.genes_exclude)
  }
  object[["genes"]] <- object[["genes"]][!bool.genes_exclude]
  print(sprintf("Number of genes in object: %s", length(object[["genes"]])))
  return(object)
}


map_to_human <- function(object, type_mouse_gene_ids) {
  gene_ids.human <- mouse_to_human_ortholog_gene_mapping(object[["genes"]], type_mouse_gene_ids) # mouse genes that did not map to human are assigned NA value.
  bool.genes_exclude <- is.na(gene_ids.human)
  for (x in c("frac_expr", "mean", "var")) {
    object[[x]] <- object[[x]] %>% filter(!bool.genes_exclude)
  }
  for (sem_name in names(object[["sem"]])) { # assumes sem and sem_bin have same elements, which is a fair assumption
    object[["sem"]][[sem_name]] <- object[["sem"]][[sem_name]] %>% filter(!bool.genes_exclude)
    object[["sem_bin"]][[sem_name]] <- object[["sem_bin"]][[sem_name]] %>% filter(!bool.genes_exclude)
  }
  object[["genes"]] <- gene_ids.human[!bool.genes_exclude] # updating genes
  print("Clearing group_by_annotation.sem and group_by_annotation.sem_bin slots. This is not important, but we do it to avoid any downstream confusion")
  object[["group_by_annotation.sem"]] <- NULL
  object[["group_by_annotation.sem_bin"]] <- NULL
  print("Info: sem, sem_bim slot are left 'as is' after mapping to human. The recommended workflow is to recalculate sem_bin slot")
  return(object)
}


calc_sem_wrapper <- function(object) {
  ### Function: wrapper to calculate all SEMs
  sems <- c("tstat", "ges", "zscore", "si", "specificity")
  for (sem in sems) {
    object <- calc_sem(object, sem_name=sem)
  }
  return(object)
}

calc_sem <- function(object, sem_name) {
  ### Function: Sets sem slot. All sems are tibbles.
  if (!is.null(object[["sem"]][[sem_name]])) { # check if slot already exists
    print(sprintf("Warning: %s slot already exists in object. Data will be overwritten", sem_name))
  }
  if (sem_name == "tstat") {
    df_sem <- tstat_sem(object)
  } else if (sem_name == "ges") {
    df_sem <- ges_sem(object)
  } else if (sem_name == "zscore") {
    df_sem <- zscore_sem(object)
  } else if (sem_name == "si") {
    df_sem <- normalized_specificity_index(object[["mean"]])
  } else if (sem_name == "specificity") {
    df_sem <- calculate_specificity(object[["mean"]])
  # } else if (sem_name == "mean") { # makes no sense to include. 
  #   df_sem <- object[["mean"]]
  # } else if (sem_name == "frac_expr") { # makes no sense to include
  #   df_sem <- object[["frac_expr"]]
  # } else if (sem_name == "pem") {
  #   df_sem <- calculate_pem(object[["mean"]])
  } else {
    stop(sprintf("Received unrecognized sem_name: %s", sem_name))
  }
  ### check for NA values
  df_sem.na <- df_sem %>% filter_all(any_vars(is.na(.)))
  if (nrow(df_sem.na) > 0) {
    print(sprintf("sem_name=%s | N=%s genes contained NA values after SEM calculation.", sem_name, nrow(df_sem.na)))
  }
  object[["sem"]][[sem_name]] <- as.tibble(df_sem)
  return(object)
}

bin_sems <- function(object, n_bins=101, threshold_bin_zero=0) {
  ### Function sets the slot 'sem_bin' containing binned data for each of the sems in the 'sem' slot
  if (is.null(object[["sem"]])) { # ensure that sem slot exists
    stop("object has no sem slot. Calculate SEMs before running this function")
  }
  if (!is.null(object[["sem_bin"]])) { # check if slot already exists
    print(sprintf("Warning: %s slot already exists in object. Data will be overwritten", "sem_bin"))
  }
  object[["sem_bin"]] <- list()
  for (name.sem in names(object[["sem"]])) {
    object[["sem_bin"]][[name.sem]] <- wrapper_binning(object[["sem"]][[name.sem]], n_bins=n_bins, threshold_bin_zero=threshold_bin_zero)
  }
  return(object)
}


set_group_by_annotation_slots <- function(object) {
  ### Function: returns a list of lists with per annotation.
  ### group sem and sim_bin data by annotations
  object[["group_by_annotation.sem"]] <- list()
  object[["group_by_annotation.sem_bin"]] <- list()
  for (annotation in object[["annotations"]]) {
    list.tmp.sem <- list()
    list.tmp.sem_bin <- list()
    for (sem_name in names(object[["sem"]])) { # assumes sem and sem_bin have same elements, which is a fair assumption
      list.tmp.sem[[sem_name]] <- object[["sem"]][[sem_name]] %>% select(!!rlang::sym(sem_name) := !!rlang::sym(annotation)) # selecting and renaming
      list.tmp.sem_bin[[sem_name]] <- object[["sem_bin"]][[sem_name]] %>% select(!!rlang::sym(sem_name) := !!rlang::sym(annotation)) # selecting and renaming
    }
    object[["group_by_annotation.sem"]][[annotation]] <- bind_cols(list.tmp.sem) # tibble
    object[["group_by_annotation.sem_bin"]][[annotation]] <- bind_cols(list.tmp.sem_bin) # tibble
  }
  return(object)
}

calc_sem_meta_from_bins <- function(object) {
  if (is.null(object[["group_by_annotation.sem_bin"]]) || (length(object[["group_by_annotation.sem_bin"]])==0) ) { # short-circuit OR
    stop("Calculate 'group_by_annotation.sem_bin' slot before running this function")
  }
  print("Calculating sem_meta based on group_by_annotation.sem_bin")
  object[["sem_meta"]] <- list()
  list.tmp.median <- list()
  list.tmp.mean <- list()
  list.tmp.sd <- list()
  for (annotation in object[["annotations"]]) {
    print(sprintf("Calculating for annotation = %s", annotation))
    list.tmp.median[[annotation]] <- apply(object[["group_by_annotation.sem_bin"]][[annotation]], 1, median)
    list.tmp.mean[[annotation]] <- apply(object[["group_by_annotation.sem_bin"]][[annotation]], 1, mean)
    list.tmp.sd[[annotation]] <- apply(object[["group_by_annotation.sem_bin"]][[annotation]], 1, sd)
  }
  object[["sem_meta"]][["median"]] <- bind_cols(list.tmp.median)
  object[["sem_meta"]][["mean"]] <- bind_cols(list.tmp.mean)
  object[["sem_meta"]][["sd"]] <- bind_cols(list.tmp.sd)
  return(object)
}


hierarchical_sem <- function(object, list.annotations) {
  ### object: 'base' object to use for hierarchical SEMs. Must contain all annotations.
  ### list.annotations: named list of annotations (character vector). The names(list.annotations) should be the names of the hierarchical grouping of annotation.
  
  ### Check that all annotations in list.annotations are in object. We do this to avoid downstream errors (eaiser then try-catch)
  for (name_elem in names(list.annotations)) {
    annotations <- list.annotations[[name_elem]]
    if (!all(annotations %in% object[["annotations"]])) {
      stop(sprintf("Error: hierarchical grouping = %s contains annotations not present in object annotations"))
    }
  }
  ### Run hierarchical loop
  list.objects <- list()
  for (name_elem in names(list.annotations)) {
    print(sprintf("Doing hierarchical SEM for hierarchical group = %s", name_elem))
    annotations <- list.annotations[[name_elem]]
    object.sub <- subset_annotations(object, annotations)
    object.sub <- calc_sem_wrapper(object.sub)
    object.sub <- bin_sems(object.sub)
    object.sub <- set_group_by_annotation_slots(object.sub)
    object.sub <- calc_sem_meta_from_bins(object.sub)
    list.objects[[name_elem]] <- object.sub
  }
  return(list.objects)
}



######################################################################################################
######################## SEM model fit #############################
######################################################################################################

library(broom)

fit_sems <- function(object, df.magma, df.metadata) {
  ### Add genetic data to sem_meta
  df.sem_meta <- object[["sem_meta"]][["median"]] %>% 
    mutate(gene=object[["genes"]]) %>%  # add genes
    inner_join(df.magma %>% select(gene, ZSTAT), by="gene") %>% # inner join
    select(gene, ZSTAT, everything())
  
  ### Run regressions
  print("Running regressions...")
  list.res <- list()
  for (annotation in object[["annotations"]]) {
    print(annotation)
    df.tmp <- df.sem_meta %>% select(ZSTAT, !!rlang::sym(annotation)) # two-column
    list.res[[annotation]] <- lm(as.formula(paste0("ZSTAT ~", annotation)), data = df.tmp)
  }
  
  df.res <- tibble::enframe(list.res, name="annotation", value="fit") # REF https://r4ds.had.co.nz/many-models.html#from-a-named-list
  
  print("Running broom...")
  df.broom <- df.res %>% mutate(
    tidied = map(fit, tidy),
    glanced = map(fit, glance)
  )
  
  df.tidied <- df.broom %>% select(-fit, -glanced) %>% unnest(tidied) %>% filter(term != "(Intercept)") # discard intercept - keep only beta for the cell-type
  df.glance <- df.broom %>% select(-fit, -tidied) %>% unnest(glanced)
  
  ### add beta and se
  df.model_sumstats <- df.glance %>% left_join(df.tidied %>% select(annotation, estimate, std.error), by="annotation") %>% arrange(p.value)
  df.model_sumstats <- df.model_sumstats %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))
  ### add metadata
  df.model_sumstats <- df.model_sumstats %>% left_join(df.metadata, by="annotation")
  
  return(df.model_sumstats)
}
  

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
######################## mouse_to_human_ortholog_gene_expression_mapping [DATA FRAME INPUT] #############################
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


