############### SYNOPSIS ###################
# Run Seurat transfer labels analysis

### OUTPUT:
#' csv table, example:
#' id_datExpr_test,cell_type_test,id_datExpr_ref,predicted.id,n_cells,median.prediction.score
#' moffitt2018.umi,e8_Cck_Ebf3,moffitt2018.umi,i9_Gaba,12,0.565629814207739

### REMARKS:
# ....

### REFERENCE:
# * 

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--dir_data", type="character", default=NULL,
              help = "character, relative path from project top dir to dir containing raw test input expression data, gene*cell table with gene names in the first column, [default %default]"), 
  make_option("--id_datExpr_test", type="character",
              help = "character, used to match raw expression matrix filename in dir_data, gene*cell table with gene names in the first column, also used as a label, e.g. 'campbell2017.umi'"), 
  make_option("--metadata_test_annot_col", type="character", default = "cell_type",
              help = "character, name of test metadata column containing cell annotation [default %default]"),
  make_option("--vec_metadata_test_annots_subset", type="character", default= NULL,
              help = "quoted vector of celltype ids to map [default $default]"), 
  make_option("--id_datExpr_ref", type="character",
              help = "character, used to match raw expression matrix filename in dir_data [default %default]"),
  make_option("--metadata_ref_annot_col", type="character", default = "cell_type",
              help = "character, name of reference metadata column containing cell annotation [default %default]"), 
  make_option("--nComp", type="integer", default= 50L,
              help = "integer, number of components to use in the reduced dimension [default $default]"), 
  make_option("--output_label", type="character",
              help = "character, label for output files"),
  make_option("--append_results", type="logical", default=TRUE,
              help = "If file already exists, append results? [default %default] ")
)


opt <- parse_args(OptionParser(option_list=option_list))

dir_data <- opt$dir_data
id_datExpr_test <- opt$id_datExpr_test
metadata_test_annot_col <- opt$metadata_test_annot_col
vec_metadata_test_annots_subset <- opt$vec_metadata_test_annots_subset
if (!is.null(opt$vec_metadata_test_annots_subset)) vec_metadata_test_annots_subset <- eval(parse(text=vec_metadata_test_annots_subset))
id_datExpr_ref <- opt$id_datExpr_ref
metadata_ref_annot_col <- opt$metadata_ref_annot_col
nComp <- opt$nComp
output_label <- opt$output_label
append_results <- opt$append_results

######################################################################
########################### PACKAGES #################################
######################################################################

message("Loading packages")

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Seurat"))

######################################################################
########################### SET OPTIONS ##############################
######################################################################

options(stringsAsFactors = F, 
        use="pairwise.complete.obs", 
        warn=1, 
        verbose=F,
        mc.cores=40 # for parallel computation
) 

######################################################################
############################ CONSTANTS ###############################
######################################################################

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

randomSeed = 12345
set.seed(randomSeed)

######################################################################
########################## LOAD SEURAT OBJECTS #######################
######################################################################

message("Loading data")

path_seuratObj_test <- dir(path=dir_data, pattern=paste0(id_datExpr_test, ".*seuratObj\\.RDS\\.gz$"), full.names=T)
seuratObj_test <- readRDS(path_seuratObj_test)
Idents(seuratObj_test) <- seuratObj_test@meta.data[[metadata_test_annot_col]]

path_seuratObj_ref <- dir(path=dir_data, pattern=paste0(id_datExpr_ref, ".*seuratObj\\.RDS\\.gz$"), full.names=T)
seuratObj_ref<- readRDS(path_seuratObj_ref)
Idents(seuratObj_ref) <- seuratObj_ref@meta.data[[metadata_ref_annot_col]]

######################################################################
###################### TRANSFER LABELS ###############################
######################################################################
# project test dataset onto reference PCA subspace

anchors <- FindTransferAnchors(reference = seuratObj_ref, 
                               reduction="pcaproject",
                               query = seuratObj_test, 
                               npcs = if (!is.null(seuratObj_ref@reductions[["pca"]])) NULL else nComp,
                               dims = 1:nComp, 
                               verbose=T)

predictions <- TransferData(anchorset = anchors, 
                            refdata = seuratObj_ref@meta.data[[metadata_ref_annot_col]],  
                            dims = 1:nComp, 
                            verbose=T) 

######################################################################
########################## prepare results ###########################
######################################################################

message("preparing outputs")

predictions <- data.table("cell_id" = rownames(predictions), predictions) 

if (!is.null(vec_metadata_test_annots_subset))  {
  seuratObj_test <- subset(x=seuratObj_test, idents=vec_metadata_test_annots_subset)
  predictions <- predictions[cell_id %in% colnames(seuratObj_test)]
}

predictions[,cell_type_test:=seuratObj_test@meta.data[[metadata_test_annot_col]]]

predictions_summary <- predictions[,.(n_cells = .N,
                                      median.prediction.score = median(prediction.score.max)),
                                   by = .(cell_type_test,predicted.id),]   
  
predictions_summary$id_datExpr_test = id_datExpr_test
predictions_summary$id_datExpr_ref = id_datExpr_ref

######################################################################
########################## prepare summary data.table ################
######################################################################

predictions_summary <- predictions_summary[order(cell_type_test, -n_cells)]

######################################################################
#############################  WRAP UP  ##############################
######################################################################

file.out.results.summary<- here(paste0(output_label, "_seurat_compare_celltypes_summary.csv"))
fwrite(x = predictions_summary, 
       file = file.out.results.summary,  
       append = if (!file.exists(file.out.results.summary)) F else append_results, 
       nThread=24, 
       verbose=T)
message(paste0("Seurat ", id_datExpr_ref, " label transfer to ",  id_datExpr_test, " done!"))


