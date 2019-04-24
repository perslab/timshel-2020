############### SYNOPSIS ###################
# Expression specificity library.


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################


library(tidyverse)


######################################################################################################
########################################## TODO IMPLEMENTATION #######################################
######################################################################################################

### general
# - consistent naming: find/replace all SEM/sem with es.
# - make into a package with roxygen documentation.

### sem_obj
# - add slot for annotation metadata (e.g. object$metadata). This will make many things simpler, because you don't have to read in metadata from separate files.

### es calculations
# - implement "exclude self" category during ES calculation

######################################################################################################
######################################## Binning ############################################
######################################################################################################

# wrapper_binning <- function(df, n_bins=101, threshold_bin_zero=0){
#   print(sprintf("Binning data frame with n_bins=%s", n_bins))
#   df.binned <- df %>% transmute_all(funs(bin_vector), n_bins, threshold_bin_zero) # *OBS* output tibble looses rownames
#   # rownames(df.binned) <- rownames(df) # set rownmaes *IMPORTANT*
#   return(df.binned)
# }

bin_vector <- function(vector_in, n_bins, threshold_bin_zero=0){
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

rank_normalize_vecor <- function(vector_in, threshold_bin_zero=-Inf){
  # INPUT: a vector of numeric values to rank normalize
  # OUTPUT: a numeric vector with values {0 ... 1}
  # *REMARK #1*: note that vector_in <= threshold_bin_zero will be given the value 0.
  # *REMARK #2*: note that NA values in vector_in will be given the value 0.
  values <- rep(0,length(vector_in))
  bool.not_bin_zero <- (vector_in > threshold_bin_zero) & (!is.na(vector_in)) # bool index for elements that should not be assigned value 0
  values[bool.not_bin_zero] <- base::rank(vector_in[bool.not_bin_zero], ties.method="average")/length(vector_in[bool.not_bin_zero])
  return(values)
}


######################################################################################################
############################# SKENE expression functions ################################
######################################################################################################


# wrapper_specificity_quantiles <- function(df.avg_expr, dataset_prefix) {
#   ### INPUT: data frame or matrix of average expression. genes x celltypes. with rownames and columnnames.
#   # SEE calculate_specificity() for notes on sensitivity to log_transformation.
#   ### OUTPUT: data frame with specificity quantiles.
#   
#   df.specificity <- calculate_specificity(df.avg_expr)
#   write_csv(df.specificity %>% rownames_to_column(var="gene"), sprintf("%s.specificity.sem.csv.gz", dataset_prefix))
#   
#   df.specificity.quantiles <- df.specificity %>% transmute_all(funs(transform_to_quantiles), threshold_bin_zero=0) # *OBS* output tibble looses rownames
#   rownames(df.specificity.quantiles) <- rownames(df.specificity) # set rownmaes *IMPORTANT*
#   write_csv(df.specificity.quantiles %>% rownames_to_column(var="gene"), sprintf("%s.specificity_quantiles.sem.csv.gz", dataset_prefix))
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
  # ^ adding 0.000000000001 avoids 'devision by zero' problems
  # ^ R divison by zero:
  # ^ 0/0=NA
  # ^ 1/0=Inf
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
############################################# SI ##############################################
######################################################################################################

# ^ XXXX then when we rank the genes while ignoring NA values (na.last="keep"), the gene will get a low rank or maybe even be NA.

normalized_specificity_index <- function (df) {
  ### Input: data.frame. genes x cell-types. Any rownames are not used
  ### Output: normalized SI data.frame with same dimensions as input df.
  ### Algorithm inspired by specificity.index() function from pSI CRAN package.
  ### pSI package references:
  # Dougherty et al. 2010: https://academic.oup.com/nar/article/38/13/4218/2409074
  # R CRAN package (v1.1, 2014-01-30, Official implementation): https://cran.r-project.org/package=pSI
  # Second implementation: http://genetics.wustl.edu/jdlab/psi_package/
  # First implementation: http://www.bactrap.org/downloads/Specificity.r
  ### Improvements of nSI (normalized specificity index) compared to the original SI
  ### Summary of differences between nSI and SI:
  # 1) Normalized ranks: if a gene is only expressed in one cell-type it will get nSI = 1. The SI score can be abitrary high.
  # 2.1) Added epsilon to ratio to avoid NA/-Inf values: relevant for single-cell data that contain many zeroes.
  # 2.2) Added numerical tricks to deal with zero values: fc=~0-->0 and fc=~1-->0
  # 2.3) Use rank(ties.method="min") to add maximal penalty to genes with fc = 0. SI uses the default rank(ties.method="average").
  # 3) We do not calculate log2(fc) because the log-transformation does not do anything [except potentially generating NA or -Inf] when ranking afterward 
  #    [NB: in the paper, the authors also do not use log2 in the SI equation (eq 1).]
  ### Additional remarks:
  # We note that the Dougherty2010 SI equation (eq.1) is NOT consistent with any of the R implementations.
  
  print("Calculating nSI")
  epsilon <- 1e-100
  df <- as.data.frame(df) # make sure the input is a data frame.
  n_annotations <- ncol(df) 
  df.SI <- data.frame(matrix(nrow = nrow(df) , ncol = n_annotations)) 
  colnames(df.SI) <- colnames(df) 
  # rownames(df.SI) <- rownames(df) # PT note: rownames are used
  for (j in 1:n_annotations) { # j: index cell-types (columns of input)
    print(sprintf("Calculating nSI: #%s/#%s", j, n_annotations))
    notme <- c(1:n_annotations)[-j] # exclude self index
    ### Calculate fold-change
    # fc <- df[, j]/df[, notme] # without epsilon
    fc <- (df[, j]+epsilon)/(df[, notme]+epsilon) # df, genes x (annotations-1) | column-wise division.
    # ^ ***MARCH 13th 2019***: I realized that we can obtain equivalent results with perhaps slightly simpler code (but less transparant) without epsilon because the way R handles division by zero 1/0=Inf and 0/0=NA
    # ^ ***MARCH 13th 2019***: Specifically, we can obtain equivalent results by removing epsilon from fc formula and not doing dplyr::near() 'mapping'.
    # ^ ***MARCH 13th 2019***: But because the code is more transparant and works on python, we keep epsilon.
    # ^ R division by zero:
      # ^ 0/0=NA
      # ^ 1/0=Inf
      # ^ -1/0=-Inf
    # ^ example to illustrate why it would work WITHOUT epsilon:
    # ^ specific_gene: umi.self=100; umi.other=0 --> fc=100/0=Inf
    # ^ unspecific_gene: umi.self=0; umi.other=10 --> fc=0/10=0
    # ^ unspecific_gene_dropouts: umi.self=0; umi.other=0 --> fc=0/0=NA
    # ^ the ranks of specific_gene and unspecific_gene will work and neeed no further explanation.
    # ^ the rank of unspecific_gene_dropouts will get a low rank IF AND ONLY IF we propagate NA values (na.last="keep") and afterwards set them to zero.
    ### With epsilon
    # ^ because we add epsilon we need to do a few tricks after calculating fc:
    # ^ fc=~0 --> 0: ok
    # ^ fc=~1 --> 0: ok, see below
    # ^ fc>>c --> do nothing: fc can get very big (e.g. ~1e100) when annotation_j > 0 and annotation_k is 0: (some_number+epsilon)/(0+epsilon). We don't need to do anything there. The ranks will still work.
    fc[dplyr::near(fc,0)] <- 0 # we map fc near 0 to 0 to avoid numerical precision artifacts when ranking the genes and keeping track of ties. fc values near 0 arise when the annotation[j] is zero and the epsilon is added: (0+epsilon)/(some_number+epsilon)
    fc[dplyr::near(fc,1)] <- 0 # we map fc near 1 to 0 because these fc values are (almost) certainly a result of dividing two zero measurements: (0+epsilon)/(0+epsilon). These fc values should be 0, because the gene is not expressed in the given annotation.
    ### Rank genes: would want high ranks to correspond to high specificity
    fc_ranks <- apply(fc, 2, base::rank, na.last="keep", ties.method="min") # matrix, (genes x annotations-1)
    # ^ rank in ascending order (low values ==> low rank; high values ==> high rank)
    # ^ na.last="keep": NA values are kept with rank NA. That is, NA values are ignored during ranking and 'propagated' to the result vector. 
    #                   we must set na.last="keep" because of how rank() treats NA values: "NA values are never considered to be equal: for na.last = TRUE and na.last = FALSE they are given distinct ranks in the order in which they occur in x."
    #                   that is, without na.last="keep" the NA values would be ARBITRARILY RANK and we don't want that when we, at least without epsilon, could get a lot of NA values in fc.
    #                   NB: na.last="keep" is really only important if we DID NOT use episilon to avoid NA values in fc.
    # ^ ties.method="min" : tied genes get the minimum rank. We do this to add maximal penalty to genes with fc = 0. 
    #                       We are awere that it will also effect genes with fc >>, but they are unlikely to have tied fc values 
    #                       emperically we observe that ties.method="min" creates a 'bimodal' distribution of nSI values: some genes will be specific, others (the majority) will be unspecific with the same value.
    # ALT FASTER: fc_ranks <- apply(fc, 2, data.table::frank, ...) # identical but faster solution using data.table::frank
    ### Normalize ranks
    fc_ranks_normalized <- rowMeans( (fc_ranks-1)/(nrow(fc_ranks)-1) ) # divide each column by the number of genes and calculate the mean.
    ### Notes on subtracting 1 in the both denominator and numerator:
    # ^ becasue we force [fc=~0-->0 and fc=~1-->0] and rank using ties.method="min", many genes will have rank=1 (lowest possible rank)
    # ^ so when we below normalize the ranks, we many genes will get the same value: fc_ranks/nrow(fc_ranks) = 1/n_genes ==> some small number.
    # ^ the fact that many genes will have nSI=1/n_genes results in a bimodal distribution of nSI [center_of_mass_no1 = 1/n_genes; center_of_mass_no2 ~ 0.8'ish],.
    # ^ By subtracting 1 in the both denominator and numerator we obtain: 
    # ^ 1) nSI will be 0 for genes with lowest possible rank. [becasue we subtract 1 in the denominator]
    # ^ 2) nSI will be 1 for genes with highest possible rank [becasue we ALSO subtract 1 in the numerator]
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

zscore_sem <- function(object_data) {
  print("Calculating z-score")
  df.zscores.gene_wise <- as.data.frame(t(scale(t(as.data.frame(object_data[["mean"]])), center = TRUE, scale = TRUE))) # returns a DATA FRAME with columns=annotations, rownames=genes.
  df.zscores.gene_wise.na <- df.zscores.gene_wise %>% filter_all(any_vars(is.na(.))) # we (sometimes) have NA values. Don't know why. It could be because some genes have zero-variance
  if (nrow(df.zscores.gene_wise.na) > 0) {
    print(sprintf("OBS: N=%s genes contained NA values after zscore calculation. ", nrow(df.zscores.gene_wise.na)))
  }
  return(df.zscores.gene_wise)
}

tstat_sem <- function(object_data, df.ncells) {
  ### vectorized calculation of pooled variance for each gene
  # formula: sum_i{var[i]*(n[i]-1)}/sum_i{n[i]-1}, where i is the group.
  
  print("Calculating t-stat")
  
  annotations <- colnames(object_data[["mean"]])
  n_genes <- dim(object_data[["mean"]])[1]
  
  df.n <- df.ncells %>% slice(rep(1:n(), each=n_genes )) # df, 1 x annotations | replicate row in ncells. REF rep(): https://stackoverflow.com/a/47780537/6639640
  sd.pooled <- sqrt(rowSums(object_data[["var"]] * (df.n-1))/rowSums(df.n - 1)) # numeric vector, genes x 1 | unbiased least squares estimate (n-1) of pooled standard deviation.
  list.tstat <- list()
  counter <- 1
  for (annotation in annotations) {
    print(sprintf("#%s/#%s | Calculating tstat for annotation=%s", counter, length(annotations), annotation))
    annotation_sym <- rlang::sym(annotation)
    x <- object_data[["mean"]] %>% pull(!!annotation_sym) # numeric vector, genes x 1 | mean self
    ### Calculate mean for 'other' using ncells as weights to recover the true mean of the 'other' group.
    # formula: sum_i{mean[i]*n[i]}/sum_i{n[i]}
    df.n_other <- df.ncells %>% select(-!!annotation_sym) %>% slice(rep(1:n(), each=n_genes )) # df, genes x annotations. Each row contains the number of cells and each row is identical.
    y <- rowSums(object_data[["mean"]] %>% select(-!!annotation_sym) * df.n_other)/rowSums(df.n) # numeric vector, genes x 1
    # ^ df.mean_other * df.n_other: column-wise multiplication. 
    # ^ rowSums(df.n %>% slice(1)) gives the same as rowSums(df.n) because all rows are identical
    n_x <- df.ncells %>% pull(!!annotation_sym) # numeric vector, genes x 1
    n_y <- rowSums(df.n_other) # numeric vector, genes x 1
    sd.pooled_scaling_factor <- sqrt(1/n_x+1/n_y) # sqrt(1/n1+1/n2)
    tstat <- (x-y)/(sd.pooled*sd.pooled_scaling_factor) # numeric, genes x 1 | tstat = (mean - other) / (pooled_sd*sqrt(1/n1+1/n2))
    # pt(tstat, df=sum(df.n)-2, lower.tail=T) # numeric vector, genes x 1 | one-sided ttest for higher expression. Each gene as the same number of degree of freedom (dof). Dof is total number of observations (cells) minus 2.
    list.tstat[[annotation]] <- tstat
    # if (annotation=="ENT1") {break}
    counter <- counter + 1
  }
  df.tstat <- bind_cols(list.tstat) # df/tbl, genes x annotation | binds by position
  return(df.tstat)
}



ges_sem <- function(object_data, df.ncells) {
  ### Calculate Gene Enrichment Score
  ### Ref: Zeisel, Cell 2018
  ### formula: (mean[celltype]+epsilon_1)/(mean[other]+epsilon_1) * (frac[celltype]+epsilon_2)/(frac[other]+epsilon_2)
  ### epsilon_1=0.1 ; epsilon_2=0.01
  print("Calculating GES")
  epsilon_1 <-  1e-100 # or .Machine$double.eps
  epsilon_2 <- 0.01
  
  annotations <- colnames(object_data[["mean"]])
  n_genes <- dim(object_data[["mean"]])[1]
  
  ### vectorized calculation of pooled variance for each gene
  # formula: sum_i{var[i]*(n[i]-1)}/sum_i{n[i]-1}, where i is the group.
  df.n <- df.ncells %>% slice(rep(1:n(), each=n_genes )) # df, 1 x annotations | replicate row in ncells. REF rep(): https://stackoverflow.com/a/47780537/6639640
  list.ges <- list()
  counter <- 1
  for (annotation in annotations) {
    print(sprintf("#%s/#%s | Calculating GES for annotation=%s", counter, length(annotations), annotation))
    annotation_sym <- rlang::sym(annotation)
    mean.self <- object_data[["mean"]] %>% pull(!!annotation_sym) # numeric vector, genes x 1
    frac.self <- object_data[["frac_expr"]] %>% pull(!!annotation_sym) # numeric vector, genes x 1
    ### Calculate mean and frac for 'other' using ncells as weights
    df.n_other <- df.ncells %>% select(-!!annotation_sym) %>% slice(rep(1:n(), each=n_genes )) # df, genes x annotations. Each row contains the number of cells and each row is identical.
    mean.other <- rowSums(object_data[["mean"]] %>% select(-!!annotation_sym) * df.n_other)/rowSums(df.n) # numeric vector, genes x 1
    frac.other <- rowSums(object_data[["frac_expr"]] %>% select(-!!annotation_sym) * df.n_other)/rowSums(df.n) # numeric vector, genes x 1
    ### Calculate GES
    ges <- (mean.self+epsilon_1)/(mean.other+epsilon_1) * (frac.self+epsilon_2)/(frac.other+epsilon_2) # numeric, genes x 1
    ges[dplyr::near(ges, 0)] <- 0 # assigning ges values 'near equal to' zero to zero.
    list.ges[[annotation]] <- ges
    counter <- counter + 1
  }
  df.ges <- bind_cols(list.ges) # df/tbl, genes x annotation | binds by position
  return(df.ges)
}


######################################################################################################
############################################# MULTI_GENESET ##############################################
######################################################################################################


write_multi_geneset_file <- function(object, 
                                     dataset_prefix, 
                                     use_raw_sem_values=F, 
                                     make_clean_annotation_names=T, 
                                     add_all_genes_in_dataset=F,
                                     es_mean_only=F, # only export object[["sem_meta"]][["mean"]]. Saves time
                                     write_file=T) {
  ### OUTPUT
  ### a file with filename multi_geneset.<dataset_prefix>.txt
  ### the file (no header) with the following columns:
  # col 1 "geneset name"  : <dataset_prefix>.<annotation=celltype>.<sem_name> or <dataset_prefix>.all_genes_in_dataset.dummy
  # col 2                       : human ensembl gene ID
  # col 3                       : value/score for the given gene. For .all_genes_in_dataset it will be all 1's.
  ### we prefix with <dataset_prefix> to allow for multiple multi_geneset files to be concatenatted together (e.g. hierarchical SEMs) and still have unique values in col1 (geneset name).
  ### only genes with sem_transformed > 0 will be written to file. This makes the LDSC make_annot.py pipeline more effective.
  if (is.null(object[["sem_transformed"]])) { # ensure that slots exists
    stop("Object has no sem_transformed slot. Calculate before running this function")
  }
  
  print(sprintf("Running with make_clean_annotation_names=%s", make_clean_annotation_names))
  
  
  if (es_mean_only & (use_raw_sem_values | add_all_genes_in_dataset)) {
    stop("es_mean_only argument cannot be used together with use_raw_sem_values or add_all_genes_in_dataset")
  }
  
  if (use_raw_sem_values) {
    print("OBS: use_raw_sem_values enabled. Will use object[['sem']] as sem value slot")
    SEM_SLOT <- "sem"
  } else {
    SEM_SLOT <- "sem_transformed"
  }
  
  list_res <- list()
  i <- 1
  if (!es_mean_only) {
    ### loop over individual sem_transformed
    for (sem_name in names(object[[SEM_SLOT]])) { 
      for (annotation in names(object[[SEM_SLOT]][[sem_name]])) {
        print(sprintf("Processing SEM = %s | counter = %s", sem_name, i))
        df.tmp <- object[[SEM_SLOT]][[sem_name]] %>% 
          select(value := !!rlang::sym(annotation)) %>% # select and rename column to a generic name so we can do bind_rows() later. LHS is the name of the column we create; RHS is the STRING of a column name in the data frame.
          mutate(.annotation = as_string(annotation), # LHS is the name of the column we create; RHS is a string
                 .sem_name = as_string(sem_name), # LHS is the name of the column we create; RHS is a string. (I KNOW THIS IS SLIGHTLY CONFUSING SYNTAX, but it works)
                 gene=object[["genes"]]) # add geneset annotation and gene names 
        ### NB: this simpler (no as_string() and .annotation) but less explicit version also works.
        # df.tmp <- object[[SEM_SLOT]][[sem_name]] %>% 
        #   select(value := !!rlang::sym(annotation)) %>% # select and rename column to a generic name so we can do bind_rows() later. LHS is the name of the column we create; RHS is the STRING of a column name in the data frame.
        #   mutate(annotation = annotation, # LHS is the name of the column we create; RHS is a string
        #          sem_name = sem_name, # LHS is the name of the column we create; RHS is a string. (I KNOW THIS IS SLIGHTLY CONFUSING SYNTAX, but it works)
        #          gene=object[["genes"]]) # add geneset annotation and gene names 
        if (!use_raw_sem_values) {
          df.tmp <- df.tmp %>% filter(value > 0) # keep only sem_transformed > 0. Do this AFTER adding gene names
        }
        list_res[[i]] <- df.tmp # df.tmp has the columns: value, geneset_name and gene.
        i <- i + 1
      }
    }
    ### add tstat_top10pct_binary [works for both conditions of use_raw_sem_values]
    ### this SEM was used in Finucane2018 (LDSC-SEG)
    sem_dependency <- "tstat"
    sem_name <- "tstat_top10pct_binary"
    if (is.null(object[[SEM_SLOT]][[sem_dependency]])) {
      stop(sprintf("object[[%s]][[%s]] is NULL. Make sure that the given SEM is calculated before running this function...", SEM_SLOT, sem_dependency))
    }
    for (annotation in names(object[[SEM_SLOT]][[sem_dependency]])) {
      print(sprintf("Processing SEM = %s | counter = %s", sem_name, i))
      df.tmp <- object[[SEM_SLOT]][[sem_dependency]] %>% 
        select(value := !!rlang::sym(annotation)) %>% # select and rename column to a generic name so we can do bind_rows() later. LHS is the name of the column we create; RHS is the STRING of a column name in the data frame.
        mutate(.annotation = as_string(annotation), # LHS is the name of the column we create; RHS is a string
               .sem_name = as_string(sem_name), # LHS is the name of the column we create; RHS is a string. (I KNOW THIS IS SLIGHTLY CONFUSING SYNTAX, but it works)
               gene=object[["genes"]]) # add geneset annotation and gene names 
      # SPECIFIC STEP
      df.tmp <- df.tmp %>% 
        filter(value > quantile(value, 0.9)) %>% # keep only genes with the TOP 10 % largest (positive) values REF: https://stackoverflow.com/a/33221186/6639640
        mutate(value = if_else(value > 0, true=1, false=0)) # binarise - this is how Skene and Finucane ran LDSC
      list_res[[i]] <- df.tmp # df.tmp has the columns: value, geneset_name and gene.
      i <- i + 1
    }
    
    ### add specificity_top10pct_binary [works for both conditions of use_raw_sem_values]
    ### this SEM was used in Skene2018
    sem_dependency <- "specificity"
    sem_name <- "specificity_top10pct_binary"
    if (is.null(object[[SEM_SLOT]][[sem_dependency]])) {
      stop(sprintf("object[[%s]][[%s]] is NULL. Make sure that the given SEM is calculated before running this function...", SEM_SLOT, sem_dependency))
    }
    for (annotation in names(object[[SEM_SLOT]][[sem_dependency]])) {
      print(sprintf("Processing SEM = %s | counter = %s", sem_name, i))
      df.tmp <- object[[SEM_SLOT]][[sem_dependency]] %>% 
        select(value := !!rlang::sym(annotation)) %>% # select and rename column to a generic name so we can do bind_rows() later. LHS is the name of the column we create; RHS is the STRING of a column name in the data frame.
        mutate(.annotation = as_string(annotation), # LHS is the name of the column we create; RHS is a string
               .sem_name = as_string(sem_name), # LHS is the name of the column we create; RHS is a string. (I KNOW THIS IS SLIGHTLY CONFUSING SYNTAX, but it works)
               gene=object[["genes"]]) # add geneset annotation and gene names 
      # SPECIFIC STEP
      df.tmp <- df.tmp %>% 
        filter(value > quantile(value, 0.9)) %>% # keep only genes with the TOP 10 % largest (positive) values REF: https://stackoverflow.com/a/33221186/6639640
        mutate(value = if_else(value > 0, true=1, false=0)) # binarise - this is how Skene and Finucane ran LDSC
      list_res[[i]] <- df.tmp # df.tmp has the columns: value, geneset_name and gene.
      i <- i + 1
    }
    ### add all genes
    if (add_all_genes_in_dataset) {
      print("Adding all genes...")
      list_res[[i]] <- tibble(value=1, 
                              gene=object[["genes"]], 
                              .annotation="all_genes_in_dataset",
                              .sem_name="dummy")
      i <- i + 1 # not needed, but good to have if you rearrange the code
    }
  } # end if (!es_mean_only)
  if (!use_raw_sem_values) {
    ### loop over sem_meta mean
    sem_name <- "sem_mean"
    for (annotation in names(object[["sem_meta"]][["mean"]])) {
      print(sprintf("Processing SEM = %s | counter = %s", sem_name, i))
      df.tmp <- object[["sem_meta"]][["mean"]] %>% 
        select(value := !!rlang::sym(annotation)) %>% # select and rename column to a generic name so we can do bind_rows() later. LHS is the name of the column we create; RHS is the STRING of a column name in the data frame.
        mutate(.annotation = as_string(annotation), # LHS is the name of the column we create; RHS is a string
               .sem_name = as_string(sem_name), # LHS is the name of the column we create; RHS is a string. (I KNOW THIS IS SLIGHTLY CONFUSING SYNTAX, but it works)
               gene=object[["genes"]]) %>% # add geneset annotation and gene names
        filter(value > 0) # keep only sem_transformed > 0. Do this AFTER adding gene names
      list_res[[i]] <- df.tmp # df.tmp has the columns: value, geneset_name and gene.
      i <- i + 1
    }
  }
  print("Doing bind_rows..")
  # rowbind
  df.multi_geneset <- bind_rows(list_res) # When row-binding, columns are matched by name
  # clean annotation name
  if (make_clean_annotation_names) {
    print("Cleaning annotation names...")
    df.multi_geneset <- df.multi_geneset %>% mutate(.annotation = clean_annotation_name(.annotation))
  }
  # make geneset_name column
  df.multi_geneset <- df.multi_geneset %>% mutate(geneset_name=paste(dataset_prefix,
                                                                     .annotation,
                                                                     .sem_name, 
                                                                     sep="."))
  # check for any NA values
  if(any(is.na(df.multi_geneset))) {
    print("Error: df.multi_geneset contains NA values. Investigate the returned data.frame for NA values, resolve the issue and rerun this function. No file will be written.")
    return(df.multi_geneset) 
  } else {
    print("Checked df.multi_geneset for NA values and found nothing. All good to go!")
  }
  
  if (write_file) {
    file.out <- sprintf("multi_geneset.%s.txt.gz", dataset_prefix)
    print(sprintf("writing file: %s", file.out))
    df.multi_geneset %>% 
      select(geneset_name, gene, value) %>% # *IMPORTANT*: re-arrange column order to agree with make_annot.py file format
      write_tsv(file.out, col_names=F) # no header
  }
  return(df.multi_geneset)
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


write_sems <- function(object, slot, dataset_prefix, dir_out) {
  ### Function: write SEM: "sem" or "sem_transformed"
  ### We recommend that dataset_prefix contains information about specie
  stopifnot(slot %in% c("sem","sem_transformed", "sem_meta", "sem_pvalues", "null"))
  print(sprintf("Writing slot=%s files for dataset_prefix=%s", slot, dataset_prefix))
  
  if (slot == "null") {
    object.data <- object[["null"]][["sem"]]
  } else {
    object.data <- object[[slot]]
  }
  
  for (name.sem in names(object.data)) {
    file.suffix <- name.sem # default value
    if (slot == "sem_transformed") {
      file.suffix <- sprintf("es_ws.%s", name.sem)
    } else if (slot == "sem_pvalues") {
      file.suffix <- sprintf("p_emp.%s", name.sem)
    } else if (slot == "sem") {
      file.suffix <- sprintf("es_w.%s", name.sem)
    } else if (slot == "null") {
      file.suffix <- sprintf("es_w_null.%s", name.sem)
    }

    file.out <- sprintf("%s/%s.%s.csv.gz", dir_out, dataset_prefix, file.suffix)
    print(sprintf("Writing file: %s", file.out))
    object.data[[name.sem]] %>% mutate(gene=object[["genes"]]) %>% select(gene, everything()) %>% write_csv(path=file.out)
    # or use data.table::fwrite() and afterwards R.utils::gzip('filename.csv',destname='filename.csv.gz')
  }
}

### OLD
# write_sem <- function(df, dataset_prefix, name.sem, ortholog_mapping=NULL) {
#   print(sprintf("Writing files for SEM = %s", name.sem))
#   write_csv(df %>% rownames_to_column(var="gene"), sprintf("%s.%s.sem.csv.gz", dataset_prefix, name.sem))
#   
#   if (!is.null(ortholog_mapping)) {
#     df.human <- mouse_to_human_ortholog_gene_expression_mapping(df, type_mouse_gene_ids=ortholog_mapping)
#     write_csv(df.human %>% rownames_to_column(var="gene"), sprintf("%s.%s.hsapiens_orthologs.sem.csv.gz", dataset_prefix, name.sem))
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
# sem_transformed --> list(sems)
# group_by_annotation.sem_transformed --> list of dfs (annotations)
# group_by_annotation.sem --> list of dfs (annotations)
# sem_meta --> list of dfs (mean,median,sd)  made from group_by_annotation.sem_transformed
# null[mean,var,frac] --> tibbles
# null[sem] --> list(sems)

clean_annotation_name <- function(annotations) {
  ### Function to clean annotation name to use for LDSC pipeline
  ### Works on both both 'scalar' and vector input
  vec.replacement <- c("/"="-", 
                       "\\s+"="_")
  return(stringr::str_replace_all(annotations, vec.replacement))
}


check_object_annotation_names <- function(annotations) {
  # annotations: character vector
  bool.white_space <- stringr::str_detect(annotations, pattern="\\s+")
  bool.fwd_slash <- stringr::str_detect(annotations, pattern="/")
  if (any(bool.white_space)) {
    print(sprintf("*WARNING*:  annotations contains whitespace (in n=%s annotations). Consider renaming annotations to improve compatability with Linux file system", sum(bool.white_space)))
  }
  if (any(bool.fwd_slash)) {
    print(sprintf("*WARNING*:  annotations contains forward slash ('/') (in n=%s annotations). Consider renaming annotations to improve compatability with Linux file system", sum(bool.fwd_slash)))
  }
  if (!any(bool.white_space,bool.fwd_slash)) {
    print("Annotation names passed check_object_annotation_names() with no remarks")
  } 
}

process_annotations_and_genes <- function(object, flag_null) {
  ### check that genes and annotations are identical
  ### **run only once because 'gene' column is removed.**
  
  ### **OBS APRIL 2019**: it is inefficient to make this data copy. It was implemented to support annotation name clean-up
  if (flag_null) {
    object_data <- object[["null"]][["data"]]
  } else {
    object_data <- object[["data"]]
  }
  
  list_annotations <- list()
  list_genes <- list()
  for (x in c("frac_expr", "mean", "var")) {
    if (anyNA(object_data[[x]])) {
      stop(sprintf("Input data contains NA values: %s", x))
    }
    list_genes[[x]] <- object_data[[x]] %>% pull(gene)
    object_data[[x]] <- object_data[[x]] %>% select(-gene) # *remove gene col*
    list_annotations[[x]] <- colnames(object_data[[x]]) # set annotations after removing gene col
  }
  if (anyNA(Reduce(f=identical_value, x=list_annotations))) {
    stop("Annotations are not identical for pre-loaded data. Input data must contain identical cell-types in the same order.")
  }
  if (anyNA(Reduce(f=identical_value, x=list_genes))) {
    stop("Genes are not identical for pre-loaded data. Input data must contain identical genes in the same order.")
  }
  
  ### since all elements are identical, we can just pick one of them
  annotations <- list_annotations[[1]]
  genes <- list_genes[[1]]
  
  ### check annotation names 
  check_object_annotation_names(annotations)
  
  if (object[["parameters"]][["clean_annotation_names"]]) {
    ### Now we know that ALL annotations match, so now we can change the names across all slots.
    # make_clean_annotation_names
    print("Making clean annotations names...")
    annotations <- clean_annotation_name(annotations)
    colnames(object_data[["frac_expr"]]) <- annotations
    colnames(object_data[["mean"]]) <- annotations
    colnames(object_data[["var"]]) <- annotations
  }

  list.res <- list("object_data"=object_data, "genes"=genes, "annotations"=annotations) 
  return(list.res)
}

preprocess_and_set_slots <- function(object) {
  # Function sets genes, annotations, data slots
  # Function check genes and annotation are identical across the different slots.
  # Gene column is removed.
  list.res <- process_annotations_and_genes(object, flag_null=FALSE)
  object[["data"]] <- list.res[["object_data"]]
  
  ### Set genes and annotations
  object[["genes"]] <- list.res[["genes"]]
  object[["annotations"]] <- list.res[["annotations"]]
  
  ### Make clean "ncells" annotation names
  if (object[["parameters"]][["clean_annotation_names"]]) {
    colnames(object[["ncells"]]) <- clean_annotation_name(colnames(object[["ncells"]]))
  }
  
  ### Last check
  stopifnot(all(object[["annotations"]] == colnames(object[["ncells"]]))) # annotation must match annotations in ncell
  
  print("Object passed all checks")
  return(object)
}

check_null <- function(object) {
  # Function checks that annotations and genes in the null is identical to the foreground.
  # Function should be called on object after preprocess_and_set_slots().
  list.res <- process_annotations_and_genes(object, flag_null=TRUE)
  object[["null"]][["data"]] <- list.res[["object_data"]]
  
  stopifnot(all(list.res[["genes"]] == object[["genes"]])) # null genes must be identical to foreground genes
  stopifnot(all(list.res[["annotations"]] == object[["annotations"]])) # null genes must be identical to foreground genes
  print("Object null passed all checks")
  return(object)
}


create_sem_object <- function(file_prefix) {
  object <- list()
  class(object) <- "sem object"
  ### Class variables/constants
  object[["parameters"]][["clean_annotation_names"]] <- TRUE
  
  ### Data
  object[["data"]][["var"]] <- read_file_fast(sprintf("%s.pre_calc.var.csv.gz", file_prefix))
  object[["data"]][["frac_expr"]] <- read_file_fast(sprintf("%s.pre_calc.frac_expr.csv.gz", file_prefix))
  object[["data"]][["mean"]] <- read_file_fast(sprintf("%s.pre_calc.mean.csv.gz", file_prefix))
  object[["genes_exclude"]] <- read_file_fast(sprintf("%s.pre_calc.sporadically_expressed_genes.anova.csv.gz", file_prefix))
  object[["ncells"]] <- read_file_fast(sprintf("%s.pre_calc.ncells.csv.gz", file_prefix))
  
  ### making 'ncells' a one row df.
  object[["ncells"]] <- object[["ncells"]] %>% 
    spread(key=annotation, value=n) %>% # 1-row tibble with annotation as colnames. spread() will re-arrange column order, so we need to fix that in the next step
    select(object[["ncells"]] %>% pull(annotation)) # df, annotations x 1 | re-arrrenge columns back to their original order.
  ### alternative WORKS.
  # n <- object[["ncells"]] %>% pull(n)
  # names(n) <- object[["ncells"]] %>% pull(annotation)
  # object[["ncells"]] <- as.tibble(as.list(n))
  
  ### Read null
  object[["null"]] <- list()
  object[["null"]][["data"]][["var"]] <- read_file_fast(sprintf("%s.pre_calc.var.null.csv.gz", file_prefix))
  object[["null"]][["data"]][["frac_expr"]] <- read_file_fast(sprintf("%s.pre_calc.frac_expr.null.csv.gz", file_prefix))
  object[["null"]][["data"]][["mean"]] <- read_file_fast(sprintf("%s.pre_calc.mean.null.csv.gz", file_prefix))
  
  ### create object slots (not needed)
  # object[["sem"]] <- list()
  # object[["sem_transformed"]] <- list()
  
  ### pre_process: this sets 
  object <- preprocess_and_set_slots(object) # function sets genes, annotations, data slots
  
  ### check null
  object <- check_null(object)
  
  return(object)
}


clear_sem_data_slots <- function(object) {
  ### Clear data slots
  print("Clearning sem data slots: sem, sem_transformed, etc...")
  object[["sem"]] <- NULL
  object[["sem_transformed"]] <- NULL
  object[["group_by_annotation.sem"]] <- NULL
  object[["group_by_annotation.sem_transformed"]] <- NULL
  object[["sem_meta"]] <- NULL
  
  object[["null"]][["sem"]] <- NULL
  return(object)
}


subset_annotation_data <- function(object_data, annotations){
  # PRIVATE FUNCTION
  list_annotations <- list()
  for (x in c("frac_expr", "mean", "var")) {
    object_data[[x]] <- object_data[[x]] %>% select(!!!rlang::syms(annotations)) # Notice !!!syms() because annotations is a vector.
    list_annotations[[x]] <- colnames(object_data[[x]])
  }
  if (anyNA(Reduce(f=identical_value, x=list_annotations))) {
    stop("Annotations are not identical for pre-loaded data. Input data must contain identical cell-types in the same order.")
  }
  return(object_data)
}


subset_annotations <- function(object, annotations) {
  ### Function: returns a object with the subset of annotations provided in 'annotations'.
  ### annotations: a character vector of annotation names to subset on. 
  ### NB: annotations are re-ordered by the order provided in annotations.
  print(sprintf("Subsetting object: keeping n=%s annotations.", length(annotations)))

  object[["data"]] <- subset_annotation_data(object[["data"]], annotations) # foreground
  object[["null"]][["data"]] <- subset_annotation_data(object[["null"]][["data"]], annotations) # null
  object[["ncells"]] <- object[["ncells"]] %>% select(!!!rlang::syms(annotations)) # ncells
  object[["annotations"]] <- object[["annotations"]][match(annotations,object[["annotations"]])] # subject annoations
  stopifnot(colnames(object[["ncells"]]) == object[["annotations"]]) # must have same annotations. We also check that the order is the same, just so everything stays in sync
  
  ### Excluding zero mean genes (after subsetting on annotations)
  bool.genes_exclude <- rowMeans(object[["data"]][["mean"]])==0
  object <- exclude_genes(object, bool.genes_exclude)
  
  ### Clear data slots
  object <- clear_sem_data_slots(object)
  return(object)
}


exclude_genes <- function(object, bool.genes_exclude){
  # Function to exclude genes
  # bool.genes_exclude: a boolean vector of genes to exclude of equal length as the number of genes in object
  
  print(sprintf("Total number of genes marked for exclusion: %s", sum(bool.genes_exclude)))
  for (x in c("frac_expr", "mean", "var")) {
    object[["data"]][[x]] <- object[["data"]][[x]] %>% filter(!bool.genes_exclude)
    object[["null"]][["data"]][[x]] <- object[["null"]][["data"]][[x]] %>% filter(!bool.genes_exclude)
  }
  object[["genes"]] <- object[["genes"]][!bool.genes_exclude]
  print(sprintf("Number of genes in object: %s", length(object[["genes"]])))
  return(object)
}

exclude_sporadic_expressed_genes <- function(object) {
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
  bool.genes_zero_mean <- rowMeans(object[["data"]][["mean"]])==0
  print(sprintf("Number of zero-mean genes: %s", sum(bool.genes_zero_mean)))
  bool.genes_exclude <- bool.genes_zero_mean | bool.sporadic_expressed 
 
  ### Exlude genes
  object <- exclude_genes(object, bool.genes_exclude)
  return(object)
}


map_to_human <- function(object, type_mouse_gene_ids) {
  gene_ids.human <- mouse_to_human_ortholog_gene_mapping(object[["genes"]], type_mouse_gene_ids) # mouse genes that did not map to human are assigned NA value.
  bool.genes_exclude <- is.na(gene_ids.human)
  object <- exclude_genes(object, bool.genes_exclude)
  object[["genes"]] <- gene_ids.human[!bool.genes_exclude] # updating genes
  
  ### Clear data slots
  object <- clear_sem_data_slots(object)
  
  ### the below code supports keeping SEM calculation data from mouse
  # for (sem_name in names(object[["sem"]])) { # assumes sem and sem_transformed have same elements, which is a fair assumption
  #   object[["sem"]][[sem_name]] <- object[["sem"]][[sem_name]] %>% filter(!bool.genes_exclude)
  #   object[["sem_transformed"]][[sem_name]] <- object[["sem_transformed"]][[sem_name]] %>% filter(!bool.genes_exclude)
  # }
  # print("Clearing group_by_annotation.sem and group_by_annotation.sem_transformed slots. This is not important, but we do it to avoid any downstream confusion")
  # object[["group_by_annotation.sem"]] <- NULL
  # object[["group_by_annotation.sem_transformed"]] <- NULL
  # print("Info: sem, sem_bim slot are left 'as is' after mapping to human. The recommended workflow is to recalculate sem_transformed slot")
  return(object)
}


calc_sem_wrapper <- function(object) {
  ### Function: wrapper to calculate all SEMs
  sems <- c("tstat", "ges", "si", "specificity") # "zscore"
  for (sem in sems) {
    object <- calc_sem(object, sem_name=sem)
    object <- calc_sem(object, sem_name=sem, null=T)
  }
  return(object)
}

calc_sem <- function(object, sem_name, null=F) {
  ### Function: Sets sem slot. All sems are tibbles.
  # if (!is.null(object[["sem"]][[sem_name]])) { # check if slot already exists
  #   print(sprintf("Warning: %s slot already exists in object. Data will be overwritten", sem_name))
  # }
  
  ### TODO: optimization inefficient to make a copy?
  if (null) {
    print("Calculating null SEM")
    object_data <- object[["null"]][["data"]]
  } else {
    object_data <- object[["data"]]
  }
  
  if (sem_name == "tstat") {
    df_sem <- tstat_sem(object_data, object[["ncells"]])
  } else if (sem_name == "ges") {
    df_sem <- ges_sem(object_data, object[["ncells"]])
  } else if (sem_name == "si") {
    df_sem <- normalized_specificity_index(object_data[["mean"]])
  } else if (sem_name == "specificity") {
    df_sem <- calculate_specificity(object_data[["mean"]])
  # } else if (sem_name == "zscore") {
  #   df_sem <- zscore_sem(object)
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
  
  if (null) {
    object[["null"]][["sem"]][[sem_name]] <- as.tibble(df_sem)
  } else {
    object[["sem"]][[sem_name]] <- as.tibble(df_sem)
  }
  return(object)
}

calc_empirical_pvalues <- function(df.obs, df.null) {
  
  ### TMP for debugging
  # object <- sem_obj.sub
  # sem_name <- "si"
  # df.obs <- object[["sem"]][[sem_name]]
  # df.null <- object[["null"]][["sem"]][[sem_name]]
  
  list.pvals <- list()
  list.qvals <- list()
  for (annotation in colnames(df.obs)) {
    ### Formula: pvals = 1-findInterval(obs, null)/nobs: finds the fraction of observations in the null that are larger than the observed values.
    # step 1: idx.insertions <- findInterval(obs, null): Find indices where elements in OBS should be inserted to maintain order in NULL.  'NULL' must be sorted (weakly) increasingly (low values --> high values)
    # step 1: idx.insertions with low indicies are OBS with low SEM values (that should be inserted in the top of NULL). idx.insertions with high indicies are OBS with high SEM values (that should be inserted in the bottom of NULL)
    # step 2: idx.insertions/nobs: finds the fraction of observations in the null that are SMALLER than the observed values
    # step 3: 1-idx.insertions/nobs: reverses the above --> finds the fraction of observations in the null that are LARGER than the observed values
    idx.insertions <- findInterval(df.obs %>% pull(!!rlang::sym(annotation)), sort(df.null %>% pull(!!rlang::sym(annotation))))
    ### findInterval REF: https://stat.ethz.ch/R-manual/R-devel/library/base/html/findInterval.html: ... This is the same computation as for the empirical distribution function, and indeed, findInterval(t, sort(X)) is identical to n * Fn(t; X[1],..,X[n]) where Fn is the empirical distribution function of X[1],..,X[n].
    # return values of findInterval(x, vec) can take on values of {0...length(vec)}
    n_null_gt_obs <- nrow(df.obs)-idx.insertions # number of null observations *greather than* our observed value. 
    # ^ e.g. if an element in idx.insertions was the most extreme in the null, we have nrow(df.obs) == idx.insertions so n_null_gt_obs=0 for that element
    pvals_empirical <- (n_null_gt_obs+1)/(nrow(df.obs)+1)
    # *OBS*: Note that the above formula includes a commonly used pseudocount (North, Curtis, and Sham, 2002; Knijnenburg et al., 2009) to avoid P-values of zero.
    ### ALTERNATIVE shorter but less transparant version: 
    # pvals_empirical <- 1-(idx.insertions/(nrow(df.obs)+1)) # *OBS*: here you should not have +1 to idx.insertions because idx.insertions==nrow(df.obs) if our observed value is the most extreme
    
    list.pvals[[annotation]] <- pvals_empirical
    list.qvals[[annotation]] <- p.adjust(pvals_empirical, method="BH")
  }
  df.pvals <- bind_cols(list.pvals) # tibble. colnames will be annotation names
  df.qvals <- bind_cols(list.qvals) # tibble. colnames will be annotation names
  list.p_and_qvals <- list("pvals" = df.pvals, "qvals" = df.qvals)
  return(list.p_and_qvals)
}


calc_empirical_pvalues_wrapper <- function(object) {
  if (is.null(object[["sem"]])) { # ensure that slots exists
    stop("Object has no sem slot. Calculate before running this function")
  }
  for (sem_name in names(object[["sem"]])) {
    print(sprintf("Calculating emperical p-values and q-values for SEM = %s", sem_name))
    list.p_and_qvals <- calc_empirical_pvalues(object[["sem"]][[sem_name]], object[["null"]][["sem"]][[sem_name]])
    object[["sem_pvalues"]][[sem_name]] <- list.p_and_qvals[["pvals"]]
    object[["sem_qvalues"]][[sem_name]] <- list.p_and_qvals[["qvals"]]
  }
  return(object)
}


get_empirical_distribution <- function(object, sem_name, annotation=NULL) {
  if (is.null(object[["sem"]][[sem_name]]) | is.null(object[["null"]][["sem"]][[sem_name]])) {
    stop(sprintf("Object has no sem_name = %s slot in either null or sem slot.", sem_name))
  }
  if (is.null(annotation)) {
    annotation <- sample(object[["annotations"]], size=1)
    print(sprintf("Chose random annotation: %s", annotation))
  } else if (!annotation %in% object[["annotations"]]) {
    stop(sprintf("Annotation %s is not in object annotations.", annotation))
  }
  # findInterval(df.obs %>% pull(!!rlang::sym(annotation)), sort(df.null %>% pull(!!rlang::sym(annotation))))
  df <- tibble(obs=object[["sem"]][[sem_name]] %>% pull(!!rlang::sym(annotation)),
                    null=object[["null"]][["sem"]][[sem_name]] %>% pull(!!rlang::sym(annotation)))
  print(summary(df))
  return(df)
}

get_empirical_pvalues_summary <- function(object, threshold, slot="sem_pvalues") {
  if (!slot %in% c("sem_pvalues", "sem_qvalues")) { # ensure that we get the correct slot argument
    stop(sprintf("Got wrong slot argument: %s", slot))
  }
  if (is.null(object[[slot]])) { # ensure that slots exists
    stop(sprintf("Object has no %s slot. Calculate them before running this function", slot))
  }
  list.summary <- list()
  for (sem_name in names(object[[slot]])) {
    list.summary[[sem_name]] <- apply(object[[slot]][[sem_name]], 2, function(col) {sum(col <= threshold)})
  }
  df.summary <- bind_cols(list.summary) # combine to tibble. dim = annotations x sems
  df.summary <- df.summary %>% mutate(annotation=colnames(object[["sem_transformed"]][[1]])) %>% select(annotation, everything()) # add annotation as column
  print(sprintf("Per-SEM summery number of genes with %s <= %s", slot, threshold))
  print(summary(df.summary))
  return(df.summary)
}

get_es_w_threshold <- function(object, threshold_pval=0.05) {
  ### Get the minimal "empirical significant" ESw for each annotation and es metric.
  
  ### DEV
  # object <- sem_obj
  # sem_name <- "si"
  # threshold_pval <- 0.05
  
  
  if (is.null(object[["sem"]]) | is.null(object[["sem_pvalues"]]) ) { # ensure that slots exists
    stop("Object has no sem or sem_pvalues slot. Calculate them before running this function")
  }
  list.res <- list()
  for (sem_name in names(object[["sem"]])) {
    df.sem <- as.data.frame(object[["sem"]][[sem_name]]) # make a data frame of ESw values so we can use [[]] indexing. *OBS*: inefficient to copy data, but damn it for now.
    df.sem[object[["sem_pvalues"]][[sem_name]] > threshold_pval] <- NA # set ESw values above threshold_pval to NA so we can later find the minimum value
    list.res[[sem_name]] <- apply(df.sem, 2, min, na.rm=T) # apply to each column
    # enframe(apply(df.sem, 2, min, na.rm=T))
    # name                                  value
    # 1 Bladder.bladder_cell                  0.761
    # 2 Bladder.bladder_urothelial_cell       0.761
    # enframe() converts named atomic vectors or lists to one- or two-column data frames
  }
  df <- bind_cols(list.res)
  # A tibble: 115 x 4
  # tstat   ges    si specificity
  # 1.88  1.76 0.761      0.0131
  # 1.93  1.81 0.761      0.0132
  df <- df %>% mutate(annotation=object[["annotations"]])
  # tstat   ges    si specificity annotation                           
  # 1  1.88  1.76 0.761      0.0131 Bladder.bladder_cell                 
  # 2  1.93  1.81 0.761      0.0132 Bladder.bladder_urothelial_cell  
  return(df)
}


transform_sems <- function(object, method, n_bins=101, threshold_pval=0.05) {
  ### Function sets the slot 'sem_transformed' containing binned data for each of the sems in the 'sem' slot
  # 
  if (is.null(object[["sem"]]) | is.null(object[["sem_pvalues"]]) ) { # ensure that slots exists
    stop("Object has no sem or sem_pvalues slot. Calculate them before running this function")
  }
  if (!is.null(object[["sem_transformed"]])) { # check if slot already exists
    print(sprintf("Warning: %s slot already exists in object. Data will be overwritten", "sem_transformed"))
  }
  accepted_method_args <- c("binning", "rank_normalize")
  if (!method %in% accepted_method_args) {
    stop(sprintf("Got wrong method argument: %s. Accepted arguments is %s", method, paste(accepted_method_args, sep=",")))
  }
  for (sem_name in names(object[["sem"]])) {
    df.sem <- as.data.frame(object[["sem"]][[sem_name]]) # copy
    df.sem[object[["sem_pvalues"]][[sem_name]] > threshold_pval] <- NA # set SEM values above threshold_pval to NA so these values go in bin zero. | this indexing works because sem_pvalues and sem have same dimensions
    if (method == "binning") {
      object[["sem_transformed"]][[sem_name]] <- as.tibble(apply(df.sem, 2, bin_vector, n_bins=n_bins))
    } else if (method == "rank_normalize") {
      object[["sem_transformed"]][[sem_name]] <- as.tibble(apply(df.sem, 2, rank_normalize_vecor))
    }
    
  }
  return(object)
}




set_group_by_annotation_slots <- function(object) {
  ### Function: returns a list of lists with per annotation.
  ### group sem and sim_bin data by annotations
  if (is.null(object[["sem"]]) | is.null(object[["sem_transformed"]]) ) { # ensure that slots exists
    stop("Object has no sem or sem_transformed slot. Calculate them before running this function")
  }
  
  object[["group_by_annotation.sem"]] <- list()
  object[["group_by_annotation.sem_transformed"]] <- list()
  for (annotation in object[["annotations"]]) {
    print(annotation)
    list.tmp.sem <- list()
    list.tmp.sem_transformed <- list()
    for (sem_name in names(object[["sem"]])) { # assumes sem and sem_transformed have same elements, which is a fair assumption
      list.tmp.sem[[sem_name]] <- object[["sem"]][[sem_name]] %>% select(!!rlang::sym(sem_name) := !!rlang::sym(annotation)) # selecting and renaming
      list.tmp.sem_transformed[[sem_name]] <- object[["sem_transformed"]][[sem_name]] %>% select(!!rlang::sym(sem_name) := !!rlang::sym(annotation)) # selecting and renaming
    }
    object[["group_by_annotation.sem"]][[annotation]] <- bind_cols(list.tmp.sem) # tibble
    object[["group_by_annotation.sem_transformed"]][[annotation]] <- bind_cols(list.tmp.sem_transformed) # tibble
  }
  return(object)
}


calc_sem_meta <- function(object) {
  if (is.null(object[["group_by_annotation.sem_transformed"]]) || (length(object[["group_by_annotation.sem_transformed"]])==0) ) { # short-circuit OR
    stop("Calculate 'group_by_annotation.sem_transformed' slot before running this function")
  }
  print("Calculating sem_meta based on group_by_annotation.sem_transformed")
  object[["sem_meta"]] <- list()
  list.tmp.median <- list()
  list.tmp.mean <- list()
  list.tmp.sd <- list()
  for (annotation in object[["annotations"]]) {
    print(sprintf("Calculating sem_meta for annotation = %s", annotation))
    list.tmp.median[[annotation]] <- apply(object[["group_by_annotation.sem_transformed"]][[annotation]], 1, median)
    list.tmp.mean[[annotation]] <- apply(object[["group_by_annotation.sem_transformed"]][[annotation]], 1, mean)
    list.tmp.sd[[annotation]] <- apply(object[["group_by_annotation.sem_transformed"]][[annotation]], 1, sd)
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
    object.sub <- calc_empirical_pvalues_wrapper(object.sub)
    object.sub <- transform_sems(object.sub, method="rank_normalize")
    object.sub <- set_group_by_annotation_slots(object.sub)
    object.sub <- calc_sem_meta(object.sub)
    list.objects[[name_elem]] <- object.sub
  }
  return(list.objects)
}



######################################################################################################
######################################### SEM model fit #############################################
######################################################################################################

library(broom)

get_sem_meta <- function(object, df.magma) {
  df.sem_meta <- object[["sem_meta"]][["median"]] %>% 
    mutate(gene=object[["genes"]]) %>%  # add genes
    inner_join(df.magma %>% select(gene, ZSTAT), by="gene") %>% # inner join
    select(gene, ZSTAT, everything())
  return(df.sem_meta)
}

get_genetic_and_sem_data <- function(object, slot, df.magma) {
  if (slot %in% names(object[["sem_meta"]])) {
    df <- object[["sem_meta"]][[slot]]
  } else if (slot %in% names(object[["sem"]])) {
    df <- object[["sem"]][[slot]]
  } else {
    names_accepted_slots <- c(names(object[["sem_meta"]]), names(object[["sem"]]))
    stop(sprintf("Got wrong slot name: %s. Accepted slots are: %s", slot, paste(names_accepted_slots, sep=", ")))
  }
  ### Add genetic data to df
  df <- df %>% 
    mutate(gene=object[["genes"]]) %>%  # add genes
    inner_join(df.magma %>% select(gene, ZSTAT), by="gene") %>% # inner join
    select(gene, ZSTAT, everything())
  return(df)
}

fit_sems <- function(object, slot, df.magma, df.metadata=NULL, exclude_bin_zero=F) {

  df.regression <- get_genetic_and_sem_data(object, slot, df.magma)
  
  ### Run regressions
  print("Running regressions...")
  list.res <- list()
  for (annotation in object[["annotations"]]) {
    print(annotation)
    df.tmp <- df.regression %>% select(ZSTAT, annotation_vals=!!rlang::sym(annotation)) # two-column. 
    # ^ We rename <annotation_col_name> to annotation_vals to ensure that the column name is syntaxically correct.
    # ^ E.g. an annotation name with a forward slash (e.g. 'a06.NG2/OPC') will be severely misinterpretated in as.formula()
    if (exclude_bin_zero) {
      # df.tmp <- df.tmp %>% filter(!!rlang::sym(annotation)!=0) # OLD DELTETE
      df.tmp <- df.tmp %>% filter(annotation_vals !=0 )
      print(sprintf("Excluding bin zero. %s genes remain in the regression...", nrow(df.tmp)))
    }
    # list.res[[annotation]] <- lm(as.formula(paste0("ZSTAT ~", annotation_vals)), data = df.tmp) # OLD DELTETE
    list.res[[annotation]] <- lm(ZSTAT ~ annotation_vals, data = df.tmp)
    # "a06.NG2/OPC"
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
  df.model_sumstats <- df.glance %>% left_join(df.tidied %>% select(annotation, estimate, std.error), by="annotation")
  
  ### correct pvals: one sided
  print("Correcting p-values")
  df.model_sumstats <- df.model_sumstats %>% mutate(p.value.orig = p.value) # save orig colum
  df.model_sumstats <- df.model_sumstats %>% mutate(p.value = p.value/2) # convert to one-sided test: summary.lm() computed a two-sided p-value; we want the one-sided
  df.model_sumstats <- df.model_sumstats %>% mutate(p.value = if_else(estimate < 0, 1-p.value, p.value)) # convert to one-sided test: find the complementary p-value for negative beta's
  df.model_sumstats <- df.model_sumstats %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))
  
  ### add metadata
  if (!is.null(df.metadata)) {
    df.model_sumstats <- df.model_sumstats %>% left_join(df.metadata, by="annotation")
  }
  ## sort by pval
  df.model_sumstats <- df.model_sumstats %>% arrange(p.value)
  return(df.model_sumstats)
}


fit_sems_tstat <- function(object, slot, df.magma, df.metadata=NULL) {
  ### Add genetic data to sem_meta
  df.regression <- get_genetic_and_sem_data(object, slot, df.magma)
  
  ### Run regressions
  print("Running regressions...")
  list.res <- list()
  for (annotation in object[["annotations"]]) {
    print(annotation)
    df.tmp <- df.regression %>% select(ZSTAT, !!rlang::sym(annotation)) # two-column
    idx.zero <- df.tmp %>% pull(!!rlang::sym(annotation)) == 0
    non_bin_zero <- df.tmp %>% filter(!idx.zero) %>% pull(ZSTAT)
    bin_zero <- df.tmp %>% filter(idx.zero) %>% pull(ZSTAT)
    list.res[[annotation]] <- t.test(non_bin_zero,bin_zero,var.equal=F,alternative="greater")
  }
  df.res <- tibble::enframe(list.res, name="annotation", value="fit") # REF https://r4ds.had.co.nz/many-models.html#from-a-named-list
  
  print("Running broom...")
  df.broom <- df.res %>% mutate(
    tidied = map(fit, tidy)
    # tidy() and glance() returns exactly the same for t.test()
  )
  df.model_sumstats <- df.broom %>% select(-fit) %>% unnest(tidied)
  df.model_sumstats <- df.model_sumstats %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))
  ### add metadata
  if (!is.null(df.metadata)) {
    df.model_sumstats <- df.model_sumstats %>% left_join(df.metadata, by="annotation")
  }
  df.model_sumstats <- df.model_sumstats %>% left_join(df.metadata, by="annotation")
  ## sort by pval
  df.model_sumstats <- df.model_sumstats %>% arrange(p.value)
  
  return(df.model_sumstats)
}

######################################################################################################
########################################### PLOTS ####################################################
######################################################################################################

boxplot_bin_zero <- function(object, slot, df.magma, df.metadata, annotation=NULL) {
  if (is.null(annotation)) {
    annotation <- sample(object[["annotations"]], size=1)
    print(sprintf("Chose random annotation: %s", annotation))
  }
  ### Add genetic data to sem_meta
  df.regression <- get_genetic_and_sem_data(object, slot, df.magma)

  df.tmp <- df.regression %>% select(ZSTAT, !!rlang::sym(annotation)) # two-column
  idx.zero <- df.tmp %>% pull(!!rlang::sym(annotation)) == 0
  df.tmp <- df.tmp %>% mutate(group = if_else(idx.zero, "bin_zero", "non_zero_bin"))
  p <- ggplot(df.tmp, aes(x=group, y=ZSTAT)) + geom_boxplot() + labs(title=annotation)
  print(p)
  return(df.tmp)
}

