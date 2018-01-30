
library(Seurat)
library(tidyverse)
library(Matrix)

library(parallel)


wd <- "/raid5/projects/timshel/sc-genetics/src/nb-avg_cell_type_expression"
setwd(wd)


### Source custom scripts
dir.project_src <- "/projects/timshel/sc-arc_lira/src"
dir.pers_lab_sc_lib <- "/projects/timshel/git/perslab-sc-library"
source(sprintf("%s/constants-sample_meta_data.R", dir.project_src)) # loads sample meta-data
source(sprintf("%s/constants-gene_markers_data.R", dir.project_src)) # loads gene markers (relevant for arc-lira)
source(sprintf("%s/constants-cell_type_annotations.R", dir.project_src)) # loads cell type annotations
source(sprintf("%s/seurat_functions/load_functions.R", dir.pers_lab_sc_lib)) # load Pers lab/Timshel single-cell library


# ======================================================================= #
# ==============================  USAGE  =============================== #
# ======================================================================= #

# time Rscript run_per_cell_type_expr_ttest.R |& run_per_cell_type_expr_ttest.out.txt

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #


# Load data
file.RData.cell_atlas <- "/projects/timshel/sc-arc_lira/src/out-generate_cell_atlas/arc_lira_cell_atlas-core_objs.RData"
load(file.RData.cell_atlas) # seurat_obj, df.cluster_markers


### global_notebook_parameters
### cell type annotation
p.res <- 4; p.res.string <- paste0("res.", p.res)
# df.cluster_annotation <- get_cell_atlas_annotation(resolution=p.res.string) # PRETTY NAMES
# df.cluster_annotation

# ## Assigning cell type identity to clusters
# seurat_obj <- SetAllIdent(seurat_obj, id=p.res.string)
# cell_type_mapped <- plyr::mapvalues(x=seurat_obj@ident, from=df.cluster_annotation$cluster, to=df.cluster_annotation$annotation)
# seurat_obj@ident <- cell_type_mapped # Set cell types ident
# seurat_obj <- StashIdent(seurat_obj, save.name = "annotation") # saves ident in @meta.data

# TSNEPlot(SetAllIdent(seurat_obj, id="annotation"), do.label=T, no.legend = TRUE, label.size=6, do.return=T)

# ======================================================================= #
# ==============================  FUNCTION  ============================== #
# ======================================================================= #

## # T-statistic cell type specific expression
ttest_specific_expression_per_cell_type <- function(colname_ident, ident, genes, seurat_obj) {
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
  t_stats <- sapply(genes, function(gene) {t.test(x = data.test[gene, cells.1], y = data.test[gene, cells.2], alternative="two.sided", var.equal=F)$statistic})
  # lapply: lapply(X, FUN, ...) - returns a list (unnamed)
  # sapply: sapply(X, FUN, ...) - returns a vector (named with gene names)
  # TODO: consider using pblapply/pbsapply for status bar
  # (Welch --> var.equal=F)
  ### Create data frame and return
  df.res <- data.frame(t_stats, row.names=genes) # add rownames
  colnames(df.res) <- ident # add column name (e.g. "Astrocyte_1")
  return(df.res)
}


# ======================================================================= #
# =================================  RUN  =============================== #
# ======================================================================= #

print("Finding markers...")
N_CORES <- 25

### Using cluster numbers
COLNAME_IDENT <- p.res.string # "res.0.8"
IDENTS <- sort(unique(as.numeric(seurat_obj@meta.data[,p.res.string]))) # numeric vector, 0, 1, 2, 3, ... (N_CLUSTERS-1)
GENES <- rownames(seurat_obj@data) # get all genes

### Using annotations
# IDENTS <- sort(df.cluster_annotation$annotation)
# COLNAME_IDENT <- "annotation"
# GENES <- rownames(seurat_obj@data) # get all genes
# # GENES <- rownames(seurat_obj@data)[1:10] # get all genes


#Identify cluster markers
cl <- makeCluster(N_CORES, type = "FORK")
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlang)) # needed for using UQ()
list_of_dfs.celltype_expr_ttest <- parLapply(cl, IDENTS, function(x_ident){ttest_specific_expression_per_cell_type(COLNAME_IDENT, x_ident, GENES, seurat_obj)})
stopCluster(cl)
df.celltype_expr_ttest <- bind_cols(list_of_dfs.celltype_expr_ttest)
# *OBS*: bind_cols(). combine list of dfs into a single data frame
# ^ bind_cols(): rows are matched *by position*, so all data frames must have the same number of rows.
rownames(df.celltype_expr_ttest) <- GENES



# ======================================================================= #
# =================================  EXPORT  =============================== #
# ======================================================================= #

file.out <- "arc_lira-celltype_expr_ttest-res4.csv.gz"
file.out.gz <- gzfile(file.out, 'w')
write.table(df.celltype_expr_ttest, file=file.out.gz, col.names=T, row.names=T, quote=F, sep=",")
close(file.out.gz)
