############### SYNOPSIS ###################
# Helper funtions for running RolyPoly

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################

library(tidyverse)
library(rlang) # needed for UQ() and sym()


######################################################################################################
############################################### GWAS DATA ############################################
######################################################################################################


read_gwas <- function(file.gwas, 
                      exlcude.HLA=T, 
                      do.log_odds=F,
                      delim="\t", # file delimiter
                      col_chrom="chr", # NB
                      col_pos="pos",
                      col_rsid="rsID", # NB
                      col_beta="beta",
                      col_se="se",
                      col_maf="snp_maf" # NB
                      ) {
  ### INPUT
  # file.gwas          file path for input GWAS. Should include a header with the above mentioned six columns.
  # do.log_odds        if the input GWAS is a case control study with odd ratios in the "beta" column, then these should be transformed to LOG ODDs, before inputting to RolyPoly.
  ### Output
  # data frame        data frame that can be directly used in RolyPoly call
  #                   with the column names: c('chrom', 'pos', 'rsid', 'beta', 'se', 'maf')
  

  print(sprintf("Reading GWAS data file: %s", file.gwas))
  df.gwas <- suppressMessages(read_delim(file.gwas, delim=delim)) # suppressMessages or use "col_types=cols()"
  
  # Rename columns to fit with RolyPoly scheme
  df.gwas.rp <- df.gwas %>% 
    rename(
      chrom=UQ(rlang::sym(col_chrom)),
      pos=UQ(rlang::sym(col_pos)),
      rsid=UQ(rlang::sym(col_rsid)),
      beta=UQ(rlang::sym(col_beta)),
      se=UQ(rlang::sym(col_se)),
      maf=UQ(rlang::sym(col_maf))
    )
  
  # Select autosomes
  df.non_autosomal_snps <- df.gwas.rp %>% filter(!chrom %in% 1:22) %>% count(chrom)
  print(sprintf("Number of non-autosomal SNPs found %s", nrow(df.non_autosomal_snps)))
  if (nrow(df.non_autosomal_snps) > 0) {
    print("Non autosomal SNPs will be excluded (example):")
    print(head(df.non_autosomal_snps))
    df.gwas.rp <- df.gwas.rp %>% filter(chrom %in% 1:22) # filtering
  }
  
  # Filter GWAS: exclude HLA
  if (exlcude.HLA) {
    ### Number of MHC/HLA SNPs [chr6:25Mb-34Mb]
    n_hla_snps <- df.gwas.rp %>% filter( (((chrom == 6)&(pos >= 25e6)) & ((chrom == 6)&(pos <= 34e6))) ) %>% nrow() # --> 13852 HLA SNPs for BMI.
    print(sprintf("Number of excluded HLA SNPs: %s", n_hla_snps))
    df.gwas.rp <- df.gwas.rp %>% filter( !(((chrom == 6)&(pos >= 25e6)) & ((chrom == 6)&(pos <= 34e6))) ) # remove MHC/HLA [chr6:25Mb-34Mb]
  }
  
  return(df.gwas.rp)
}




######################################################################################################
########################################### EXPRESSION DATA #########################################
######################################################################################################
# TODO: change name pos_transformation to numerical_transformation?

read_expression_data <- function(file.expr,
                                 scale_genes, 
                                 pos_transformation="square", 
                                 genes.filter=NULL,  
                                 col_genes="gene", 
                                 delim=","
                                 ) {
  ### INPUT
  # file.expr               file path for annotation specific expression data. Columns = cells/tissues. Rows = genes.
  #                         file must contain (human ensembl) genes in the column specified by "col_geness"
  # scale_genes             boolean. If true, the gene values will be *SCALED* (zscores).
  #                         we recommend using this option when 'average expression profiles' are used to quantify annotation specific expression.
  # pos_transformation      numerical transformation used to make the specific expression values POSITIVE (semi-required by RolyPoly). OPTIONS:
  #                         "square": squared values (default)
  #                         "abs": absolute values
  #                         "pos_only": replace negative values with zero.
  #                         "none": do not transform any values.
  # genes.filter            vector with gene IDs to use for filtering. E.g. genes from 'df.gene_annot'
  # col_genes               column name for the gene IDs.
  ### Output
  # data frame              data frame that can be directly used in RolyPoly call
  
  # ROLYPOLY REQUIREMENTS: Each column name is a cell-type/tissue names. Rownames must be genes (and match the 'label' from the gene annotation)
  
  print(sprintf("Reading gene expression file: %s", file.expr))
  df.expr <- suppressMessages(read_delim(file.expr, delim=delim))
  
  if (!col_genes %in% names(df.expr)) {
    print(names(df.expr))
    stop(sprintf("col_genes=%s was not found in the expression data"))
  }
  
  ### Filter
  if (!is.null(genes.filter)) {
    bool.expr_genes_keep <- df.expr %>% pull(UQ(rlang::sym(col_genes))) %in% genes.filter
    print("Genes filtering enabled.")
    print(sprintf("Number of genes filtered out: %s (genes not contained in genes.filter)", nrow(df.expr) - sum(bool.expr_genes_keep)))
    print(sprintf("Number of genes kept: %s", sum(bool.expr_genes_keep)))
    df.expr <- df.expr %>% filter(bool.expr_genes_keep)
  }
  
  ## Format for RolyPoly: set "col_genes" as rownames.

  df.expr <- df.expr %>% 
    as.data.frame() %>%
    column_to_rownames(var = col_genes)
    #column_to_rownames(var = UQ(rlang::sym(col_genes)) ) # OBS: here the UQ(rlang:sym()) does NOT work. This is likely because we have called "as.data.frame()"
  

  if (scale_genes) {
    print("Will scale data across genes (mean=0 and sd=1 for each gene)..")
    df.expr.t <- t(as.data.frame(df.expr)) # transposing. returns a DATA FRAME with columns=genes, rownames=annotations
    df.expr.t.scale <- scale(df.expr.t, center = TRUE, scale = TRUE) # returns a *MATRIX* with colnumns=genes and rownames=annotations
    df.expr <- as.data.frame(t(df.expr.t.scale)) # transposing back. returns a DATA FRAME with columns=annotations, rownames=genes.
    #ALTERNATIVE (less memory usage):
    # df.expr <- as.data.frame(t(scale(df.expr.t, center = TRUE, scale = TRUE)))
  } else {
    print("Will not scale genes.")
  }
  
  
  ### Numerical transformation of expression data: convert to positive values
  print(sprintf("Transforming annotation specific expression values using transformation: %s", pos_transformation))
  if (pos_transformation=="square") {
    df.expr.rp <- df.expr^2 # *OBS*: x^2 or x**2 returns a matrix!
  } else if (pos_transformation=="abs") {
    df.expr.rp <- abs(df.expr) # returns data.frame
  } else if (pos_transformation=="pos_only") {
    df.expr.rp <- replace(df.expr, df.expr < 0, 0) # replace negative values with 0
  } else if (pos_transformation=="none") {
    df.expr.rp <- df.expr # no transformation
  } else {
    stop(sprintf("Error: got wrong 'pos_transformation' argument: %s", pos_transformation))
  }
  
  df.expr.rp <- as.data.frame(df.expr.rp) # making sure we return a data frame
  return(df.expr.rp)
}


######################################################################################################
############################################ ANNOTATIONS DATA ########################################
######################################################################################################

### SNIPPET file.gene_annot
# ensembl_gene_id chromosome_name start_position end_position   gene_biotype
# 1 ENSG00000261657     HG991_PATCH       66119285     66465398 protein_coding
# 2 ENSG00000223116              13       23551994     23552136          miRNA
# 3 ENSG00000233440              13       23708313     23708703     pseudogene
# 4 ENSG00000207157              13       23726725     23726825       misc_RNA

read_gene_annotation <- function(file.gene_annot, protein_coding_only=F, delim="\t") {
  print(sprintf("Reading gene annotation file: %s", file.gene_annot))
  df.gene_annot <- suppressMessages(read_delim(file.gene_annot, delim=delim)) # suppressMessages or use "col_types=cols()"
  
  ### Filter genes: filter out non-autosomal genes (discard X, Y, contigs etc)
  n.discarded <- df.gene_annot %>% filter(!chromosome_name %in% 1:22) %>% nrow()
  print(sprintf("Number of non-autosomal genes discarded: %s", n.discarded))
  df.gene_annot <- df.gene_annot %>% filter(chromosome_name %in% 1:22)
  
  ### Filter genes: protein coding
  if (protein_coding_only) {
    print("Option 'protein_coding_only' enabled: keeping only protein_coding genes")
    df.gene_annot <- df.gene_annot %>% filter(gene_biotype == "protein_coding")
  }
  print(sprintf("Final number of genes in df.gene_annot: %s", nrow(df.gene_annot)))
  
  return(df.gene_annot)
}


