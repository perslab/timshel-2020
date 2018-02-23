############### SYNOPSIS ###################
# Analyzing multiple RolyPoly run

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################



# ======================================================================================================= #
#=============================================== USAGE ================================================== #
# ======================================================================================================= #



# ======================================================================================================= #
# ============================================ OptParse ================================================= #
# ======================================================================================================= #

library(optparse)
option_list <- list( 
  make_option("--input_file", type="character", default=NULL,
              help = "File path to RData object file")
)
opt <- parse_args(OptionParser(option_list=option_list))

file.rp <- opt$input_file

# file.rp <- "out.rolypoly_objs.body_BMI_Locke2015.squared_tss_10kb.final.RData"
#file.rp <- "out.rolypoly_objs-v1/out.rolypoly_objs.disease_CARDIOVASCULAR.squared_tss_10kb.final.RData"

name.gwas <- stringr::str_split(basename(file.rp), "\\.")[[1]][3] # "out.rolypoly_objs.body_BMI_Locke2015.squared_tss_10kb.final.RData" --->  name.gwas = "body_BMI_Locke2015"
print(name.gwas)

# ======================================================================================================= #
# ==========================================       CONFIG      =========================================== #
# ======================================================================================================= #

# wd <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/"
# setwd(wd)

# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #

library(rolypoly)
library(tidyverse)
library(ggplot2)

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# names.expr_data.incl_tstat <- c("arc_lira.avg_expr",
#   "arc_lira.ttest",
#   "gtex.sub_tissue.avg_expr",
#   "maca.per_celltype.avg_expr",
#   "maca.per_celltype.ttest",
#   "maca.per_tissue_celltype.avg_expr",
#   "maca.per_tissue_celltype.ttest",
#   "maca.per_tissue.avg_expr",
#   "maca.per_tissue.ttest")


names.expr_data.incl_tstat <- c("arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
  "arc_lira/arc_lira.celltype_expr.ttest.hsapiens_orthologs.csv.gz",
  "gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
  "maca/maca.per_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
  "maca/maca.per_celltype.celltype_expr.ttest.hsapiens_orthologs.csv.gz",
  "maca/maca.per_tissue_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
  "maca/maca.per_tissue_celltype.celltype_expr.ttest.hsapiens_orthologs.csv.gz",
  "maca/maca.per_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
  "maca/maca.per_tissue.celltype_expr.ttest.hsapiens_orthologs.csv.gz")

# names.expr_data.avg_expr <- c("arc_lira.avg_expr",
#   "gtex.sub_tissue.avg_expr",
#   "gtex.tissue.gene_tpm.avg_expr",
#   "maca.per_celltype.avg_expr",
#   "maca.per_tissue_celltype.avg_expr",
#   "maca.per_tissue.avg_expr")


# names.expr_data.avg_expr <- c("arc_lira/arc_lira.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
#                               "gtex/gtex.sub_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
#                               "gtex/gtex.tissue.gene_tpm.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
#                               "maca/maca.per_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
#                               "maca/maca.per_tissue_celltype.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz",
#                               "maca/maca.per_tissue.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz")

# ======================================================================================================= #
# =========================================== LOAD ROLY POLY ============================================ #
# ======================================================================================================= #


# file.rp <- "out.rolypoly_objs.TEST_body_BMI_Locke2015.gene.10kb.squared.protein_coding_only.inference.RData"; name.gwas <- "TEST_body_BMI_Locke2015"
# file.rp <- "out.rolypoly_objs.TEST_body_BMI_Locke2015.gene.10kb.squared.protein_coding_only.gwas_linked.RData"; name.gwas <- "TEST_body_BMI_Locke2015"

# file.rp <- "out.rolypoly_objs.disease_CARDIOVASCULAR.squared_tss_10kb.final.RData"
# file.rp <- "out.rolypoly_objs-v1/out.rolypoly_objs.body_BMI_Locke2015.squared_tss_10kb.final.RData"; name.gwas <- "body_BMI_Locke2015"
# file.rp <- "out.rolypoly_objs.lipids_TC.Willer2013.squared_tss_10kb.final.RData"; name.gwas <- "lipids_TC_Willer2013"
# file.rp <- "out.rolypoly_objs.cov_EDU_YEARS.squared_tss_10kb.final.RData"; name.gwas <- "cov_EDU_YEARS"
# file.rp <- "out.rolypoly_objs.body_WHRadjBMIz.squared_tss_10kb.final.RData"; name.gwas <- "body_WHRadjBMIz"
# file.rp <- "out.rolypoly_objs.disease_T2D.squared_tss_10kb.final.RData"; name.gwas <- "disease_T2D"
# file.rp <- "out.rolypoly_objs.body_HEIGHTz.squared_tss_10kb.final.RData"; name.gwas <- "body_HEIGHTz"
# file.rp <- "out.rolypoly_objs.mental_SCZ_Ripke2014.squared_tss_10kb.final.RData"; name.gwas <- "mental_SCZ_Ripke2014"

tmp.loadtime <- system.time(load(file.rp)) # rp.first_run, runtime.first_run, list.rp_wrapper
print(tmp.loadtime)
names.list.rp_wrapper.orig <- names(list.rp_wrapper) # save orig names

if (length(list.rp_wrapper) == 9) { # incl. tstat
  names(list.rp_wrapper) <- names.expr_data.incl_tstat
} else if (length(list.rp_wrapper) == 6) { # avg_expr
  # names(list.rp_wrapper) <- names.expr_data.avg_expr
  # do nothing
} else {
  stop("Wrong list.rp_wrapper length")
}


# ======================================================================================================= #
# ============================================= CLEAN OBJECT ============================================ #
# ======================================================================================================= #

GWAS_NAME <- name.gwas
RUN_NAME <- "squared_tss_10kb"

file.out.gwas_linked <- sprintf("out.rolypoly_objs-v1.clean/out.rolypoly_objs.%s.%s.gwas_linked.RData", GWAS_NAME, RUN_NAME)
file.out.inference <- sprintf("out.rolypoly_objs-v1.clean/out.rolypoly_objs.%s.%s.inference.RData", GWAS_NAME, RUN_NAME)


### Clean list.rp_wrapper
for (name.expr.data in names(list.rp_wrapper)) {
  list.rp_wrapper[[name.expr.data]]$rp$data <- "REMOVED"
  # list.rp_wrapper[[name.expr.data]]$rp$blocks <- "REMOVED"
  # list.rp_wrapper[[name.expr.data]]$rp$raw_block_data <- "REMOVED"
}
list.rp_inference <- list.rp_wrapper
### Save
save(list.rp_inference, list.run_parameters, file=file.out.inference)

### Save gwas_lined
list.gwas_linked <- list("rp"=list(rp.first_run), "runtime"=list(runtime.first_run))
save(list.gwas_linked, list.run_parameters, file=file.out.gwas_linked)



# ======================================================================================================= #
# ============================================ SIMPLE STATS ============================================= #
# ======================================================================================================= #

### new data structure
# lapply(lapply(list.rp_wrapper, '[[', 'runtime'), '[[', 'elapsed')


# ======================================================================================================= #
# -============================================== SET RP =============================================== #
# ======================================================================================================= #


# rp <- list.rp_wrapper[["gtex.sub_tissue.avg_expr"]]$rp
# rp <- list.rp_wrapper[["maca.per_celltype.ttest"]]$rp
# rp <- list.rp_wrapper[["maca.per_celltype.avg_expr"]]$rp
# 

### Things to extract
# rp$block_heritability_contribution (h_k^annot?)


# for (name.expr in names(list.rp_wrapper)) {
#   dir.out <- sprintf("export.%s", name.gwas) # no trailing slash
#   if (dir.out) { # remove any existing folder
#     unlink(dir.out)
#   } else {
#     dir.create(dir.out,showWarnings=F)
#   }
# 
#   print(sprintf("Processing %s", name.expr))
#   
#   ### Get rp object
#   rp <- list.rp_wrapper[[name.expr]]$rp
#   
#   ### Plot the log10(p) to rank tissues by the strength of their association
#   p.pval_ranking <- plot_rolypoly_annotation_ranking(rp) + labs(title=name.expr) + theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=50.5, unit="points"))
#   file.out.plot <- sprintf("%s/%s.plot_pval_ranking.%s.pdf", dir.out, name.gwas, name.expr)
#   ggsave(p.pval_ranking, filename=file.out.plot, width=30 , height=10)
#   
#   ### Plot the TOP 10 log10(p) to rank tissues by the strength of their association
#   
#   top10_annotations <- rp$bootstrap_results %>% 
#     filter(annotation != "intercept") %>% 
#     arrange(-bt_value) %>% 
#     slice(1:10) %>% 
#     pull(annotation) %>% 
#     as.character()
#   p.pval_ranking <- plot_rolypoly_annotation_ranking(rp) + labs(title=name.expr) + theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=50.5, unit="points"))
#   p.pval_ranking <- p.pval_ranking + xlim(top10_annotations)
#   file.out.plot <- sprintf("%s/%s.plot_pval_top10_ranking.%s.pdf", dir.out, name.gwas, name.expr)
#   ggsave(p.pval_ranking, filename=file.out.plot, width=14, height=7)
#   
#   ### Plot gamma estimate with 95% confidence intervals
#   p.gamma_est <- plot_rolypoly_annotation_estimates(rp) + labs(title=name.expr, x="gamma_hat_k") + theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=50.5, unit="points"))
#   file.out.plot <- sprintf("%s/%s.plot_gamma_estimates.%s.pdf", dir.out, name.gwas, name.expr)
#   ggsave(p.gamma_est, filename=file.out.plot, width=30, height=10)
#   
#   ### Export P-values
#   file.out.table <- sprintf("%s/%s.table_pval_ranking.%s.csv", dir.out, name.gwas, name.expr)
#   write_csv(rp$bootstrap_results %>% arrange(-bt_value), path=file.out.table)
#   
#   ### Export gamma-estimates
#   df.gamma <- rp$full_results$parameters %>% # "rp$full_results$parameters" is a named character vector
#     as_tibble() %>% # as_tibble/as.tibble/as_data_frame are aliases.
#     rownames_to_column(var = "annotation") %>% # for some reason, you must call rownames_to_column() before calling arrange() - else you loose the rownames.
#     arrange(desc(value))
#   file.out.table <- sprintf("%s/%s.table_gamma.%s.csv", dir.out, name.gwas, name.expr)
#   write_csv(df.gamma, path=file.out.table)
# }



# ======================================================================================================= #
# -============================================== PLOTS =============================================== #
# ======================================================================================================= #
# 
# rp <- rp.first_run
# 
# ### Plot the log10(p) to rank tissues by the strength of their association
# p.pval_ranking <- plot_rolypoly_annotation_ranking(rp)
# p.pval_ranking
# 
# ### Plot gamma estimate with 95% confidence intervals
# ### gamma: influence of cell-type-specific gene expression on the variance of GWAS effect sizes
# p.gamma_est <- plot_rolypoly_annotation_estimates(rp)
# p.gamma_est
# 
# 
# rp$full_results$parameters
# 
# names(rp$block_heritability_contribution)


# ======================================================================================================= #
# ========================================= Calculate h_j_gene ========================================== #
# ======================================================================================================= #



# rp.rolypoly_gene_score <- rolypoly_add_ld_corrected_gwas_block_scores(rolypoly=rp, fast_calculation=T) # this takes time
       # sets the slot rp$data$<BLOCK/GENE>$corrected_gwas_block_score_pval
### You might get the warning:
# In CompQuadForm::imhof(raw_score, lambda = lambda, epsabs = 1e-15,  ... :
#                          Note that Qq + abserr is positive.

# df.h_j_gene <- sapply(rp.rolypoly_gene_score$data, '[[', "corrected_gwas_block_score_pval") %>% # extract h_j_gene | returns named vector with p-values
#   as_tibble() %>% # as_tibble/as.tibble/as_data_frame are aliases.
#   rownames_to_column(var = "h_j_gene") %>% # for some reason, you must call rownames_to_column() before calling arrange() - else you loose the rownames.
#   arrange(value)



# ======================================================================================================= #
# ======================================== RP OBJECT SIZE =============================================== #
# ======================================================================================================= #

### Method for calculating memory used by objects
# rp <- list.rp_wrapper[["maca.per_tissue.avg_expr"]]$rp # maca.per_tissue.avg_expr / maca.per_tissue_celltype.avg_expr
# res <- lapply(rp, object.size)
# res
# res %>% as_tibble() %>% t()/1000000
  # blocks                                      5.048280
  # data                                      605.841400
  # raw_block_data                              3.188640
  # full_results                                0.001856
  # block_values                                1.024248
  # block_heritability_contribution             0.001520
  # expected_block_values                       1.111176
  # bootstrap_results                           0.005760
  # bootstrap_block_values                      2.675920
  # bootstrap_block_heritability_contribution   0.005576


# ======================================================================================================= #
# ======================================== RP OBJECT STRUCTURE ========================================== #
# ======================================================================================================= #

### Method for clearing data slots
# rp$data <- "REMOVED" # gene level information (GWAS-linked) [takes up 98% of the size of the object]
# rp$blocks <- "REMOVED" # gene==block annotation
# rp$raw_block_data <- "REMOVED" # gene==block expression data

### rp$data$ENSG00000000457 --> gene level information (GWAS-linked)
  # block_info ---> gene annotation (single row)
  # snps
  # ld_data --> data frame with columns: SNP_A, SNP_B, R
  # n_snps --> number of SNPs
  # gwas_block_score --> sum(beta_squared)
  # ld_matrix_squared --> [n_snps x n_snps] sparse Matrix of class "dsCMatrix"
  # partition
  # y
  # x
  # include_in_inference

### rp$blocks --> block==gene annotation
# label chrom     start       end partition
# 1     ENSG00000223972     1      1869     21869         1
# 2     ENSG00000227232     1      4363     24363         1
# 3     ENSG00000243485     1     19554     39554         1

### rp$raw_block_data --> gene==block expression data
# Aorta      Bladder Brain_Microglia Brain_Non.microglia  
# ENSG00000081791 4.335426e-03 3.715975e-01    9.470096e-01        2.542704e+00 
# ENSG00000162929 7.021306e-02 5.935347e-01    9.480770e-01        1.717095e+00 
# ENSG00000168887 5.211482e-01 1.323443e+00    1.652628e+00        7.242094e-04 


# ======================================================================================================= #
# ================================================ MISC ================================================= #
# ======================================================================================================= #


# rp$full_results$parameters %>% sort
# rp$bootstrap_results %>% arrange(-bt_value) %>% head



# ======================================================================================================= #
# ============================================== XXXXXXXX =============================================== #
# ======================================================================================================= #


