
# =========================================================================== #
# ============================= SUMMARY OF FIXES ============================ #
# =========================================================================== #

### format_sumstats_for_magma.r + standardise.sumstats.column.headers.r
# 1. Changed sed commands: removed empty qoutes. FROM "sed -i '' '1s/%s/%s/'" to "sed -i '1s/%s/%s/' %s"
# 2. Changed gawk commands. "gawk -i inplace" requires GNU Awk >= v4.1.0 (2013 release). Redhat only has GNU Awk 4.0.2

### map.snps.to.genes.r
# [Semi hacky] Fixed 'bug' to make magma --annotate command use .bim file in /tools/magma/1.06/g1000_eur.
# [Semi hacky --> deleted] Changed magma --pval command to use --bfile '%s/g1000_eur'

### generate.celltype.data.r
# Added no_cores argument

# =========================================================================== #
#========================= generate.celltype.data.r ========================= #
# =========================================================================== #

#' generate.celltype.data
#'
#' \code{generate.celltype.data} Takes expression & cell type annotations and creates celltype_data files which contain the mean and specificity matrices
#'
#' @param exp Numerical matrix with row for each gene and column for each cell. Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be cross referenced against the annot data frame.
#' @param annotLevels List with arrays of strings containing the cell type names associated with each column in exp
#' @param groupName A human readable name for refering to the dataset being loaded
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("cortex_mrna")
#' expData = cortex_mrna$exp
#' l1=cortex_mrna$annot$level1class
#' l2=cortex_mrna$annot$level2class
#' annotLevels = list(l1=l1,l2=l2)
#' fNames_ALLCELLS = generate.celltype.data(exp=expData,annotLevels,"allKImouse")
#' @export
#' @import parallel
generate.celltype.data.PT <- function(exp,annotLevels,groupName,no_cores){
  require("parallel")
  # Calculate the number of cores
  # no_cores <- detectCores() # PT outcommented
  cl <- makeCluster(no_cores)
  
  # First, check the number of annotations equals the number of columns in the expression data
  lapply(annotLevels,test <- function(x,exp){if(length(x)!=dim(exp)[2]){stop("Error: length of all annotation levels must equal the number of columns in exp matrix")}},exp)
  
  # Check group name
  if(is.null(groupName)){stop("ERROR: groupName must be set. groupName is used to label the files created by this function.")}
  if(groupName==""){stop("ERROR: groupName must be set. groupName is used to label the files created by this function.")}
  
  # Convert characters to numbers
  exp2<-suppressWarnings(apply(exp,2,function(x) {storage.mode(x) <- 'double'; x}))
  
  ctd = list()
  for(i in 1:length(annotLevels)){ctd[[length(ctd)+1]] = list(annot=annotLevels[[i]])}
  
  aggregate.over.celltypes <- function(rowOfMeans,celltypes){
    exp_out = as.matrix(data.frame(aggregate(rowOfMeans,by=list(celltypes),FUN=mean)))
    rownames(exp_out) = exp_out[,"Group.1"]
    exp_out = exp_out[,2]
    exp_out2 = as.numeric(exp_out)
    names(exp_out2) = names(exp_out)
    return(exp_out2)
  }
  calculate.meanexp.for.level <- function(ctd_oneLevel,expMatrix){
    if(dim(expMatrix)[2]==length(unique(ctd_oneLevel$annot))){
      print(dim(expMatrix)[2])
      print(length(ctd_oneLevel$annot))
      if(sum(!colnames(expMatrix)==ctd_oneLevel$annot)!=0){
        stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
      }
      ctd_oneLevel$mean_exp = expMatrix
    }else{
      mean_exp = apply(expMatrix,1,aggregate.over.celltypes,ctd_oneLevel$annot)
      ctd_oneLevel$mean_exp = t(mean_exp)
    }
    return(ctd_oneLevel)
  }
  calculate.specificity.for.level <- function(ctd_oneLevel){
    normalised_meanExp = t(t(ctd_oneLevel$mean_exp)*(1/colSums(ctd_oneLevel$mean_exp)))
    ctd_oneLevel$specificity = normalised_meanExp/(apply(normalised_meanExp,1,sum)+0.000000000001)
    return(ctd_oneLevel)
  }
  print("Running step 1: mclapply(ctd,calculate.meanexp.for.level,exp2)...")
  ctd2 = mclapply(ctd,calculate.meanexp.for.level,exp2)
  print("Running step 2: mclapply(ctd2,calculate.specificity.for.level)...")
  ctd3 = mclapply(ctd2,calculate.specificity.for.level)
  ctd=ctd3
  stopCluster(cl)
  
  fNames=sprintf("CellTypeData_%s.rda",groupName)
  print("Done. Now saving...")
  save(ctd,file=fNames)
  return(fNames)
}

# =========================================================================== #
# ================= map.snps.to.genes.r =================== #
# =========================================================================== #

#' Map SNPs to their nearby genes
#'
#' Make two external calls to MAGMA. First use it to annotate SNPs onto their neighbouring genes. Second, use it to calculate the gene level trait association.
#'
#' @param path Filepath of the summary statistics file
#' @param upstream_kb How many kb upstream of the gene should SNPs be included?
#' @param downstream_kb How many kb downstream of the gene should SNPs be included?
#' @param N What is the N number for this GWAS? That is cases+controls
#' @param genome_ref_path Path to the folder containing the 1000 genomes .bed files (which can be downloaded from https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip)
#'
#' @return Filepath for the genes.out file
#'
#' @examples
#' genesOutPath = map.snps.to.genes(gwas_sumstats_path)
#'
#' @export
map.snps.to.genes.PT <- function(gwas_sumstats_path,upstream_kb=10,downstream_kb=1.5,N=NULL,genome_ref_path){
  gwas_sumstats_path = path.expand(gwas_sumstats_path)
  
  # Check whether there is an N column in the sumstats file (if it wasn't provided as an argument)
  if(is.null(N)){
    con <- file(gwas_sumstats_path,"r") ; first_line <- readLines(con,n=1) ; close(con)
    column_headers = strsplit(first_line,"\t")[[1]]
    if("N" %in% column_headers){n_arg = "ncol=N"}else{
      nval <- as.numeric(readline("There is no N column within the sumstats file. What is the N value for this GWAS?"))
      
      if(is.na(nval)){stop(sprintf("%s provided but value of N for the GWAS must be numeric",nval))}
      if(nval<1000){stop("Value of N provided is less than 1000. This seems unlikely.")}
      if(nval>100000000){stop("Value of N provided is over than 100000000. In 2018 this seems unlikely.")}
      n_arg = sprintf("N=%s",nval)
    }
  }
  
  # Determine which genome build it uses & get path to gene loc file
  genome_build = get_genomebuild_for_sumstats(gwas_sumstats_path)
  gene_loc_dir = sprintf("%s/data/",system.file(package="MAGMA.Celltyping"))
  if(genome_build == "GRCh37"){genomeLocFile=sprintf("%s/NCBI37.3.gene.loc",gene_loc_dir)}
  if(genome_build == "GRCh38"){genomeLocFile=sprintf("%s/NCBI38.gene.loc",gene_loc_dir)}
  print(sprintf("GWAS Sumstats appear to come from genome build: %s",genome_build))
  
  sumstatsPrefix = sprintf("%s.%sUP.%sDOWN",gwas_sumstats_path,upstream_kb,downstream_kb)
  # magma_cmd = sprintf("magma --annotate window=%s,%s --snp-loc '%s' --gene-loc '%s' --out '%s'",upstream_kb,downstream_kb,gwas_sumstats_path,genomeLocFile,sumstatsPrefix)
  magma_cmd = sprintf("magma --annotate window=%s,%s --snp-loc '%s.bim' --gene-loc '%s' --out '%s'",upstream_kb,downstream_kb,path.expand(genome_ref_path),genomeLocFile,sumstatsPrefix) # *PT Edit*
  # ^ PT : The SNP location file should contain three columns: SNP ID, chromosome, and base pair position. These should be the first three columns in that file (additional columns are ignored). 
  # ^ PT : The only exception is if you use a .bim file from binary PLINK data, this can be provided to MAGMA as a SNP location file without modification.
  system(magma_cmd)
  
  # SCHIZ CLOZUK N=35802
  # magma_cmd = sprintf("magma --bfile '%s/' --pval '%s' %s --gene-annot '%s.genes.annot' --out '%s'",path.expand(genome_ref_path),gwas_sumstats_path,n_arg,sumstatsPrefix,sumstatsPrefix)
  magma_cmd = sprintf("magma --bfile '%s' --pval '%s' %s --gene-annot '%s.genes.annot' --out '%s'",path.expand(genome_ref_path),gwas_sumstats_path,n_arg,sumstatsPrefix,sumstatsPrefix) # *PT EDIT*
  # ^ PT: This will perform gene analysis on GWAS data in binary PLINK format ([DATA].bed/.bim/.fam files) using the previously generated annotation file. 
  magma_cmd
  system(magma_cmd)
  
  # Return path to genes.out file
  return(sprintf("%s.genes.out",gwas_sumstats_path))
}

# =========================================================================== #
# ================= standardise.sumstats.column.headers.r =================== #
# =========================================================================== #

#' Standardise the column headers in the Summary Statistics files
#'
#' Use a reference data table of common column header names (stored in sumstatsColHeaders.rda) convert them to a standard set, i.e. chromosome --> CHR
#' 
#' This function does not check that all the required column headers are present
#' 
#' The amended header is written directly back into the file
#'
#' @param path Filepath for the summary statistics file
#'
#' @return The amended column headers (also the column headers will be written directly into the summary statistics file)
#'
#' @examples
#' col_headers = standardise.sumstats.column.headers("~/Downloads/202040.assoc.tsv")
#'
#' @export
standardise.sumstats.column.headers.PT <- function(path){    
  # Check the sumstats file exists
  if(!file.exists(path)){stop("Path to GWAS sumstats is not valid")}
  
  # Read in the first line of the file only
  con <- file(path,"r") ; first_line <- readLines(con,n=1) ; close(con)
  column_headers = strsplit(first_line,"\t")[[1]]
  
  # Show to the user what the header is
  print("First line of summary statistics file: ")
  print(first_line)
  print(column_headers)
  
  # Amend the column headers based on a data table of commonly used names
  data(sumstatsColHeaders)
  column_headers = toupper(column_headers)
  for(headerI in 1:dim(sumstatsColHeaders)[1]){
    un = sumstatsColHeaders[headerI,1]
    cr = sumstatsColHeaders[headerI,2]
    #print(un)
    if(un %in% column_headers & (!cr %in% column_headers)){column_headers=gsub(sprintf("^%s$",un),cr,column_headers)}       
    #if(tolower(un) %in% column_headers){column_headers=gsub(sprintf("^%s$",tolower(un)),cr,column_headers)}               
  }
  new_first_line = paste(column_headers,collapse = "\t")
  
  # Write the new column headers to file
  # sed_command = sprintf("sed -i '' '1s/%s/%s/' %s",first_line,new_first_line,path) # *ORIG*
  sed_command = sprintf("sed -i '1s/%s/%s/' %s",first_line,new_first_line,path) # *PT EDIT*
  system2("/bin/bash", args = c("-c", shQuote(sed_command)))
  
  #column_headers = strsplit(column_headers," ")[[1]]
  return(column_headers)
}



# =========================================================================== #
# ======================= format_sumstats_for_magma.r ======================= #
# =========================================================================== #

#' Check that sumstats has correct columns and that they are in the correct order for MAGMA and LDSC
#'
#' @return col_headers The new column headers for the sumstats file
#'
#' @examples
#' format_sumstats_for_magma(path)
#'
#' @import data.table
#' @import stringr
#' @export
format_sumstats_for_magma.PT <- function(path){
  # Check the sumstats file exists
  if(!file.exists(path)){stop("Path to GWAS sumstats is not valid")}
  
  # col_headers = standardise.sumstats.column.headers(path)
  col_headers = standardise.sumstats.column.headers.PT(path) # *PT EDIT*
  
  
  # Check if there are CHR and BP columns
  if(!sum(c("SNP","BP") %in% col_headers)==2){
    # If not, see if there is a column storing the data in a format like: CHR:BP:A2:A1 or just CHR:BP
    # - UKBB data from Ben Neale has a Variant column with CHR:POS:REF:ALT where ALT allele is the effect allele in the model [NB: the ALT allele is NOT always the minor allele]
    # -- For input to LDSC, A1 is effect allele, A2 is non-effect allele
    # - DIAGRAM diabetes data has a Chr:Position column with CHR:BP
    # - BMI adjusted for smoking has markername with CHR:BP (with the chromosome name having 'chr' preceeding)
    # - Agression [EAGLE] just doesn't have any CHR or BP data
    print("Summary statistics file does not have obvious CHR or BP columns. Checking to see if they are joined in another column")
    
    # Obtain a row of the actual data
    con <- file(path,"r") ; row_of_data <- strsplit(readLines(con,n=2)[2],"\t")[[1]] ; close(con)
    
    # Check if there is a column of data with CHR:BP:A2:A1 format
    fourStepCol = grep(".*:.*:\\w:\\w",row_of_data)
    if(length(fourStepCol)){
      # Convert the ':' into '\t'
      # awkSplitCmd = sprintf("gawk -i inplace -F\":\" '$1=$1' OFS=\"\t\" %s",path,path)
      awkSplitCmd = sprintf("gawk -F\":\" '$1=$1' OFS=\"\t\" %s > tmp_gawk_file && mv tmp_gawk_file %s",path,path,path) # *PT EDIT*
      system2("/bin/bash", args = c("-c", shQuote(awkSplitCmd)))
      # Replace the column name with four names
      curColName = col_headers[fourStepCol]
      # Write the new column headers to file
      first_line = paste(col_headers,collapse = "\t")
      new_first_line = gsub(curColName,"CHR\tBP\tA2\tA1",paste(col_headers,collapse = "\t"))
      #sed_command = sprintf("sed -i '' '1s/%s/%s/' %s",first_line,new_first_line,path)
      sed_command = sprintf("sed -i '1s/%s/%s/' %s",first_line,new_first_line,path) # *PT EDIT*
      system2("/bin/bash", args = c("-c", shQuote(sed_command)))
      col_headers = strsplit(new_first_line,"\t")[[1]]
      print(sprintf("Column %s has been replaced with CHR BP A2 A1",curColName))
      print(col_headers)
      con <- file(path,"r") ; row_of_data <- strsplit(readLines(con,n=2)[2],"\t")[[1]] ; close(con)
    }
    
    # Check if there is a column of data with CHR:BP format
    twoStepCol = grep(".*:.*",row_of_data)
    if(length(twoStepCol)){
      # Convert the ':' into '\t'
      # awkSplitCmd = sprintf("gawk -i inplace -F\":\" '$1=$1' OFS=\"\t\" %s",path,path)
      awkSplitCmd = sprintf("gawk -F\":\" '$1=$1' OFS=\"\t\" %s > tmp_gawk_file && mv tmp_gawk_file %s",path,path,path) # *PT EDIT*
      system2("/bin/bash", args = c("-c", shQuote(awkSplitCmd)))
      # Replace the column name with four names
      curColName = col_headers[twoStepCol]
      # Write the new column headers to file
      first_line = paste(col_headers,collapse = "\t")
      new_first_line = gsub(curColName,"CHR\tBP",paste(col_headers,collapse = "\t"))
      # sed_command = sprintf("sed -i '' '1s/%s/%s/' %s",first_line,new_first_line,path)
      sed_command = sprintf("sed -i '1s/%s/%s/' %s",first_line,new_first_line,path) # *PT EDIT*
      system2("/bin/bash", args = c("-c", shQuote(sed_command)))
      col_headers = strsplit(new_first_line,"\t")[[1]]
      print(sprintf("Column %s has been replaced with CHR BP",curColName))
      print(col_headers)
      con <- file(path,"r") ; row_of_data <- strsplit(readLines(con,n=2)[2],"\t")[[1]] ; close(con)
    }
    
    # Restandardise incase the joined column headers were unusual
    col_headers = standardise.sumstats.column.headers(path)
  }
  
  # If CHR and BP are present... BUT not SNP then need to find the relevant SNP ids
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con); col_headers = strsplit(rows_of_data[1],"\t")[[1]]
  if(sum(c("CHR","BP") %in% col_headers)==2 & sum("SNP" %in% col_headers)==0){
    print("There is no SNP column found within the data. It must be inferred from CHR and BP information.")
    print("Note: this process drops any SNPs which are not from Hapmap")
    genomebuild <- as.numeric(readline("Which genome build is the data from? 1 for GRCh37, 2 for GRCh38"))
    if(!genomebuild %in% c(1,2)){stop("Genome build must be entered as either 1 (for GRCh37) or 2 (for GRCh38)")}
    data("SNP_LOC_DATA")
    if(genomebuild==1){genomebuild="GRCh37"}else{genomebuild="GRCh38"}
    snpLocDat = SNP_LOC_DATA[SNP_LOC_DATA$Build==genomebuild,][,-4]
    library(data.table)
    sumstats = fread(path)
    sumstats$CHR = as.factor(sumstats$CHR)
    if(length(grep("chr",sumstats$CHR[1]))!=0){sumstats$CHR = gsub("chr","",sumstats$CHR)}
    sumstats2 = merge(sumstats,snpLocDat,by=c("CHR","BP"))
    fwrite(sumstats2,file=path,sep="\t")
  }
  
  # Check that all the vital columns are present
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con); col_headers = strsplit(rows_of_data[1],"\t")[[1]]
  for(key_column in c("SNP","CHR","BP","P","A1","A2")){
    # code_example = "sed -i '' '1s/p_value/P/' IQ.Sniekers.2017.txt"
    code_example = "sed -i '1s/p_value/P/' IQ.Sniekers.2017.txt" # *PT EDIT*
    if(!key_column %in% col_headers){
      print("Header of file:")
      #system(sprintf("head -n 3 %s",path))
      print(rows_of_data)
      stop(sprintf("cannot find a %s column in GWAS sumstats file. \nUse code such as '%s' to fix",key_column,code_example))
    }
  }
  
  # Check there is at least one signed sumstats column
  signed_stat_column_names = c("Z","OR","BETA","LOG_ODDS","SIGNED_SUMSTAT")
  if(sum(signed_stat_column_names %in% col_headers)<1 %in% col_headers){
    print("Header of file:")
    print(rows_of_data)
    stop("ERROR: cannot find a column name representing signed statistic in GWAS sumstats file. I.e. Z, OR, BETA")
  }
  
  # Check that first three column headers are SNP, CHR, BP (in that order)
  if(!sum(col_headers[1:3]==c("SNP","CHR","BP"))==3){
    whichSNP = which(col_headers=="SNP")[1]
    whichCHR = which(col_headers=="CHR")[1]
    whichBP = which(col_headers=="BP")[1]
    otherCols = setdiff(1:length(col_headers),c(whichSNP,whichCHR,whichBP))
    #system(sprintf("gawk -i inplace '{ print $%s \" \" $%s \" \" $%s}' %s",whichSNP,whichCHR,whichBP,path))
    #newColOrder = sprintf("$%s",paste(c(whichSNP,whichCHR,whichBP,otherCols),collapse = " \" \" $"))
    newColOrder = sprintf("$%s",paste(c(whichSNP,whichCHR,whichBP,otherCols),collapse = " \"\\t\" $"))
    
    # system(sprintf("gawk -i inplace '{ print %s}' %s",newColOrder,path))
    system(sprintf("gawk '{ print %s}' %s > tmp_gawk_file && mv tmp_gawk_file %s",newColOrder,path,path)) # *PT EDIT*
  }
  
  # The formatting process can (rarely) result in duplicated columns, i.e. CHR, if CHR:BP is expanded and one already exists... delete duplicates
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con); col_headers = strsplit(rows_of_data[1],"\t")[[1]]
  if(sum(duplicated(col_headers))>0){
    notDup = which(!duplicated(col_headers))
    newColOrder = sprintf("$%s",paste(notDup,collapse = " \"\\t\" $"))
    #system(sprintf("gawk -i inplace '{ print %s}' %s",newColOrder,path))
    system(sprintf("gawk '{ print %s}' %s > tmp_gawk_file && mv tmp_gawk_file %s",newColOrder,path,path)) # *PT EDIT*
  }
  
  # MAGMA cannot handle P-values as low as 3e-400... so convert them to zeros
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con); col_headers = strsplit(rows_of_data[1],"\t")[[1]]
  # shCmd = sprintf("gawk -i inplace '{ i=%s; if($i > 1) { $i=0; }  print }' %s",which(col_headers=="P"),path.expand(path))
  shCmd = sprintf("gawk '{ i=%s; if($i > 1) { $i=0; }  print }' %s > tmp_gawk_file && mv tmp_gawk_file %s",which(col_headers=="P"),path.expand(path), path.expand(path)) # PT EDIT
  system2("/bin/bash", args = c("-c", shQuote(shCmd)))
  # The above command converts tabs in the header to spaces... so revert that (if neccesary)
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con)
  #if(("2" %in% grep("\t",rows_of_data)) & length(grep("\t",rows_of_data))==1){
  if(length(grep("\t",rows_of_data))!=2){
    tmpName = tempfile()
    shCmd = sprintf("tr ' ' '\t' < %s <> %s > %s",path.expand(path),path.expand(path),tmpName)
    system2("/bin/bash", args = c("-c", shQuote(shCmd)))
    shCmd = sprintf("mv %s %s",tmpName,path.expand(path))
    system2("/bin/bash", args = c("-c", shQuote(shCmd)))
  }
  # The above command converts 'P' to '0'... so revert that
  # shCmd = sprintf("sed -i '' '1s/0/P/' '%s'",path.expand(path))
  shCmd = sprintf("sed -i '1s/0/P/' '%s'",path.expand(path)) # *PT EDIT*
  system2("/bin/bash", args = c("-c", shQuote(shCmd)))
  
  
  # Sometimes the N column is not all integers... so round it up
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con); col_headers = strsplit(rows_of_data[1],"\t")[[1]]
  if("N" %in% col_headers){
    whichN = which(col_headers %in% "N")
    #shCmd = sprintf("gawk -i inplace '{OFS=FS="\t"}NR>1{$%s=sprintf("%3.0f",$%s)}1' %s",whichN,whichN,path.expand(path))
    # pt1 = "gawk -i inplace '{OFS=FS=\"\t\"}NR>1{$" # *ORIG*
    pt1 = "gawk '{OFS=FS=\"\t\"}NR>1{$" # *PT EDIT*
    pt2 = whichN
    pt3 = "=sprintf(\"%3.0f\",$"
    pt4 = whichN
    pt5 = ")}1' "
    # pt6 = path.expand(path) # *ORIG*
    pt6 = sprintf("%s > tmp_gawk_file && mv tmp_gawk_file %s", path.expand(path), path.expand(path)) # *PT EDIT*
    shCmd = paste(c(pt1,pt2,pt3,pt4,pt5,pt6),collapse="") # *ORIG*
    system2("/bin/bash", args = c("-c", shQuote(shCmd)))
  }
  
  # All rows should start with either SNP or rs... if they don't drop them
  shCmd = sprintf("grep -e '^rs' -e '^SNP' '%s' > '%s_tmp'",path.expand(path),path.expand(path))
  system2("/bin/bash", args = c("-c", shQuote(shCmd)))
  shCmd = sprintf("mv %s_tmp %s",path.expand(path),path.expand(path))
  system2("/bin/bash", args = c("-c", shQuote(shCmd)))     
  
  # Keep only rows which have the number of columns expected
  # shCmd = sprintf("gawk -i inplace -F'\t' 'NF == %s {print}' '%s'",length(col_headers),path.expand(path)) # *ORIG*
  shCmd = sprintf("gawk -F'\t' 'NF == %s {print}' '%s' > tmp_gawk_file && mv tmp_gawk_file %s",length(col_headers),path.expand(path), path.expand(path)) # *PT EDIT*
  system2("/bin/bash", args = c("-c", shQuote(shCmd)))
  
  # Show how the data now looks
  print("Succesfully finished preparing sumstats file:")
  print("Header of file:")
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con)
  print(rows_of_data)
  col_headers = strsplit(rows_of_data[1],"\t")[[1]]
  
  # Drop any rows with duplicated RSIDs
  # shCmd = sprintf("gawk -i inplace '!seen[$1]++' '%s'",path.expand(path)) # *ORGI*
  shCmd = sprintf("gawk '!seen[$1]++' '%s' > tmp_gawk_file && mv tmp_gawk_file %s",path.expand(path), path.expand(path)) # *PT EDIT*
  system2("/bin/bash", args = c("-c", shQuote(shCmd)))
  
  return(col_headers)
}


#all_sumstats <- readLines(path,n=-1)
#first_line = all_sumstats[1]
