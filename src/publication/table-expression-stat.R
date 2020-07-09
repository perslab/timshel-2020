############### SYNOPSIS ###################
# produce annotation * gene table with statistic of gene expression for some annotation

### OUTPUT:
#' csv table, example:

### REMARKS:
# ....

### REFERENCE:
# * 

######################################################################
########################### PACKAGES #################################
######################################################################

library("magrittr")
library("here")
library("data.table")

######################################################################
########################### OptParse #################################
######################################################################

# HARD CODE OPTIONS
opt=list(
  dir_data = here("tmp-data","expression"),
  regex_datExpr = "campbell2017\\.umi\\.csv\\.gz",
  dir_metadata = here("tmp-data","expression"),
  regex_metadata = "campbell2017\\.metadata\\.csv",
  metadata_cell_id_col = "cell_id",
  metadata_annot_col = "cell_type_all_lvl2",
  feature_subset  = fread(here("data", "genes_obesity","prot_mono_early_extr_obesity_23genes.csv"))[["gene_symbol"]],
  stat = "proportion_expr",
  output_label = "table-campbell_obesitygenes_prop_expr",
  dir_out = here("out","publication", "tables")
)



if (F) {
  suppressPackageStartupMessages(library(optparse))
  
  option_list <- list(
    make_option("--dir_data", type="character", 
                help = "character, absolute path to dir containing raw test input expression data, [default %default]"), 
    make_option("--regex_datExpr", type="character",
                help = "character, used to match raw expression matrix filename in dir_data, gene*cell table with gene names in the first column, also used as a label, e.g. 'campbell2017.umi'"), 
    make_option("--dir_metadata", type="character", 
                help = "character, absolute path to dir containing raw test input expression data, [default %default]"), 
    make_option("--regex_metadata", type="character",
                help = "character, used to match raw expression matrix filename in dir_data, gene*cell table with gene names in the first column, also used as a label, e.g. 'campbell2017.umi'"), 
    make_option("--metadata_cell_id_col", type="character", default = "cell_id",
                help = "character, name of test metadata column containing cell ids [default %default]"),
    make_option("--metadata_annot_col", type="character", default = "cell_type",
                help = "character, name of test metadata column containing cell annotation [default %default]"),
    make_option("--feature_subset", type="character", default= NULL,
                help = "whether to output statistic for a subset of features. If provided, provide either a path to a RDS file, or a quoted vector: ''c('feature1','feature2',..)'',[default $default]"), 
    make_option("--stat", type="character", default= "percentile",
                help = "output value - one of 'percentile', 'proportion_expr' [default $default]"), 
    make_option("--output_label", type="character",
                help = "character, label for output files"),
    make_option("--dir_out", type="character",
                help = "character, where to output files")
    
  )
  
  opt <- parse_args(OptionParser(option_list=option_list))
}

dir_data <- opt$dir_data
regex_datExpr <- opt$regex_datExpr
dir_metadata <- opt$dir_metadata
regex_metadata <- opt$regex_metadata
metadata_cell_id_col <- opt$metadata_cell_id_col
metadata_annot_col <- opt$metadata_annot_col
stat <- opt$stat
feature_subset <- opt$feature_subset
if (any(grepl('\\"', feature_subset))|any(grepl("\\'", feature_subset))) {
  feature_subset <- eval(parse(text=feature_subset))
} 
output_label <- opt$output_label
dir_out <- opt$dir_out



######################################################################
########################### SOURCE FNCS ##############################
######################################################################

#source(here("perslab-sc-library","functions_sc.R"))

######################################################################
########################### LOAD DATA ################################
######################################################################

path_data = dir(path = dir_data, pattern = regex_datExpr, full.names = T)

if (length(path_data)>1) {
  stop("regex_datExpr matches more than one file in dir_data")
} else if (length(path_data)==0) {
  stop("regex_datExpr matches no files in dir_data")
}


path_metadata = dir(path = dir_metadata, pattern = regex_metadata, full.names = T)
if (length(path_metadata)>1) stop("regex_metadata matches more than one file in dir_metadata")

if (!is.null(feature_subset)) {
  if (any(grepl("/",feature_subset))) {
    feature_subset = readRDS(feature_subset)
  } 
}

dt_data = fread(path_data)

dt_metadata = fread(path_metadata)

if (ncol(dt_data)-1 != nrow(dt_metadata)) stop("ncol(dt_data) + 1 must be equal to nrow(dt_metadata)")
cell_order = match(colnames(dt_data)[-1], dt_metadata[[metadata_cell_id_col]])

if (any(is.na(cell_order))) stop("dt_metadata metadata_cell_id_col entries do not all match columns of dt_data")

if (!all(feature_subset %in% dt_data[[1]])) {
  feats_missing = feature_subset[!feature_subset %in% dt_data[[1]]]
  warning(paste0(paste0(feats_missing, collapse= ", "), " not found in data"))
  feature_subset = feature_subset[!feature_subset %in% feats_missing]
}

dt_metadata = dt_metadata[cell_order,]

df_metadata = as.data.frame(dt_metadata)

rownames(df_metadata) = dt_metadata[[metadata_cell_id_col]]

vec_annot = dt_metadata[[metadata_annot_col]]

######################################################################
######################## COMPUTE OUTPUT ##############################
######################################################################

if (stat=="percentile") {
  mat_data = as.matrix(dt_data[,2:ncol(dt_data)])
  rownames(mat_data) = dt_data[[1]]
  seuratObj = CreateSeuratObject(counts = mat_data, meta.data = df_metadata)
  seuratObj = NormalizeData(seuratObj, normalization.method = "LogNormalize")
  Idents(seuratObj) = seuratObj[[metadata_annot_col]]
  vec_neuron_ids = unique(grep("^n",Idents(seuratObj),value=T))
  seuratObj = subset(x=seuratObj, idents = vec_neuron_ids)
  seuratObj_avg = Seurat::AverageExpression(object=seuratObj, slot="data", return.seurat = T)
  
  compute_percentiles = function(matIn) {
    # given a gene*sample matrix, compute percentiles
    mat_out=apply(MARGIN=2, X=matIn, FUN = function(eachcol) {
      feat_order = order(as.numeric(eachcol))
      feat_percentiles = feat_order/length(eachcol)
      names(feat_percentiles) = rownames(matIn)
      return(feat_percentiles)
    })
    
    return(mat_out)
    
  }

  if (F) {
    seuratObj_avg_counts = Seurat::AverageExpression(object=seuratObj, slot="counts", return.seurat = T)
    
    dt_data_lognorm = data.table("cell_id"=colnames(seuratObj),"annotation"= seuratObj[[metadata_annot_col]], t(as.matrix(seuratObj@assays$RNA@data)))
    colnames(dt_data_lognorm)[2] = "annotation"
    dt_data_lognorm_avg = dt_data_lognorm[,lapply(.SD,function(eachcol) {log(mean(expm1(eachcol))+1)}),by=annotation,.SDcols=colnames(dt_data_lognorm)[3:ncol(dt_data_lognorm)], ]
    mat_data_lognorm_avg = as.matrix(dt_data_lognorm_avg[,-1])
    rownames(mat_data_lognorm_avg) = dt_data_lognorm_avg$annotation
    mat_out_lognorm_manual = compute_percentiles(matIn = t(mat_data_lognorm_avg))
    mat_out_lognorm_manual = mat_out_lognorm_manual[feature_subset,]
    mat_out_lognorm_manual = mat_out_lognorm_manual[,sort(colnames(mat_out_lognorm_manual))]
    
    dt_data_counts = data.table("cell_id"=colnames(seuratObj),"annotation"= seuratObj[[metadata_annot_col]], t(as.matrix(seuratObj@assays$RNA@counts)))
    colnames(dt_data_counts)[2] = "annotation"
    dt_data_counts_avg = dt_data_counts[,lapply(.SD,function(eachcol) {log(mean(eachcol)+1)}),by=annotation,.SDcols=colnames(dt_data_counts)[3:ncol(dt_data_counts)], ]
    mat_data_counts_avg = as.matrix(dt_data_counts_avg[,-1])
    rownames(mat_data_counts_avg) = dt_data_counts_avg$annotation
    mat_out_counts_manual = compute_percentiles(matIn = t(mat_data_counts_avg))
    mat_out_counts_manual = mat_out_counts_manual[feature_subset,]
    mat_out_counts_manual = mat_out_counts_manual[,sort(colnames(mat_out_counts_manual))]
    
    mat_out_counts_seurat = compute_percentiles(matIn = as.matrix(seuratObj_avg_counts@assays$RNA@counts))
    mat_out_counts_seurat = mat_out_counts_seurat[feature_subset,]
    mat_out_counts_seurat = mat_out_counts_seurat[,sort(colnames(mat_out_counts_seurat))]
  }
  mat_out = compute_percentiles(matIn = as.matrix(seuratObj_avg@assays$RNA@data))
  mat_out = mat_out[feature_subset,]
  mat_out = t(mat_out[,sort(colnames(mat_out))])
  
} else if (stat=="proportion_expr") {
  dt_data_sub = dt_data[dt_data[[1]] %in% feature_subset]
  mat_data_sub = as.matrix(dt_data_sub[,2:ncol(dt_data_sub)])
  rownames(mat_data_sub) = dt_data_sub[[1]]
  mat_data_sub = mat_data_sub[feature_subset,]
  mat_out = apply(X = mat_data_sub, MARGIN=1, FUN=function(geneRow) {
    sapply(sort(unique(vec_annot)), function(annot) {
      sum(as.numeric(geneRow)[vec_annot==annot]>0)/sum(vec_annot==annot)
    })
  }) 
  
}

######################################################################
######################## WRITE OUT RESULTS ###########################
######################################################################



flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

dt_out = data.table("annotation"=rownames(mat_out), mat_out)
fwrite(x=dt_out, file = paste0(dir_out, "/", output_label, ".csv"))
#openxlsx::write.xlsx(x=dt_out, file = paste0(dir_out, output_label, "_", flag_date, ".xlsx"))

message("script done!")
