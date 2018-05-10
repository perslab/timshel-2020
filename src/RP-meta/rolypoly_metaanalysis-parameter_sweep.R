############### SYNOPSIS ###################
# Analyzing RolyPoly results across multiple RolyPoly parameters

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
# ==========================================       CONFIG      =========================================== #
# ======================================================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/"
setwd(wd)

# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #

suppressMessages(library(rolypoly))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))

suppressMessages(library(GGally)) # for ggpairs and ggcorr


dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================================================= #
# =========================================== LOAD ROLYPOLY ============================================ #
# ======================================================================================================= #


### Get list of files
# dir.rolypoly_objs <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-parameter_sweep/" 
# dir.rolypoly_objs <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v1.clean/" # GWAS
dir.rolypoly_objs <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v2.pos_only/" # GWAS

list.rolypoly_obj_files <- file.path(dir.rolypoly_objs, list.files(dir.rolypoly_objs, pattern="*.inference.RData")) # absolute path
head(list.rolypoly_obj_files)
print(sprintf("Got %s files to load", length(list.rolypoly_obj_files)))


### LOOP
list.gwas_inference <- list() # initialyzing empty list
for (file.rp in list.rolypoly_obj_files) {
  
  ### GWAS
  # file.rp <- list.rolypoly_obj_files[[1]] # debug
  name.gwas <- stringr::str_split(basename(file.rp), "\\.")[[1]][2] 
    # "out.rolypoly_objs.body_BMI_Locke2015.squared_tss_10kb.final.RData" --->  name.gwas = "body_BMI_Locke2015"
  
  ### PARAMETER SWEEP
  # name.gwas <- stringr::str_match(basename(file.rp), "out.rolypoly_objs.(.*).inference.RData")[[2]]
    # ^ out.rolypoly_objs.body_BMI_Locke2015.gene.10kb.squared.all_genes.inference.RData --> name.gwas == "body_BMI_Locke2015.gene.10kb.squared.all_genes"
  
  print(sprintf("loading GWAS = %s ...", name.gwas))
  tmp.loadtime <- system.time(load(file.rp)) # list.rp_inference, list.run_parameters
  print(tmp.loadtime)
  
  list.gwas_inference[[name.gwas]] <- list("inference"=list.rp_inference, "parameters"=list.run_parameters, "gwas"=name.gwas)
}


# list.gwas_inference[[<GWAS>]][["inference"]][[<EXPR_DATASET>]][["rp"]]


# ======================================================================================================= #
# =============================================== TMP SAVE ============================================== #
# ======================================================================================================= #

### GWAS
# tmp.time <- system.time(save(list.gwas_inference, file="data-list.gwas_inference-GWAS.RData"))
# print("*************************************************** SAVE TIME ***************************************************")
# print(tmp.time)
# load(file="data-list.gwas_inference-GWAS.RData")

### PARAM SWEEP
# save(list.gwas_inference, file="data-list.gwas_inference-param_sweep.RData")
# load(file="data-list.gwas_inference-param_sweep.RData")

# ======================================================================================================= #
# ====================================== SET NAMES (expression data) ==================================== #
# ======================================================================================================= #

for (name.gwas in names(list.gwas_inference)) {
  names(list.gwas_inference[[name.gwas]][["inference"]]) <- str_replace_all(names(list.gwas_inference[[name.gwas]][["inference"]]), c("^.*/(.*?)\\." = "\\1.", 
                                                                                                                                        "celltype_expr." = "", 
                                                                                                                                        ".hsapiens_orthologs.csv.gz" = ""))
}
names(list.gwas_inference[[1]][["inference"]])

# ======================================================================================================= #
# -============================================== TMP RP =============================================== #
# ======================================================================================================= #


# rp <- list.gwas_inference[["body_BMI_Locke2015.gene.10kb.squared.protein_coding_only"]][["inference"]][["gtex.tissue.gene_tpm.avg_expr"]]$rp
# df.1 <- rp$bootstrap_results
# rp$full_results$parameters


# ======================================================================================================= #
# ======================================= EXTRACT AND CREATE DATA FRAME ================================= #
# ======================================================================================================= #

### Extract results
list.dfs_meta_expr <- list()
for (name.gwas in names(list.gwas_inference)) {
  for (name.expr in names(list.gwas_inference[[name.gwas]][["inference"]])) {
    ### Create a one row data frame
    tmp.annotations <- list.gwas_inference[[name.gwas]][["inference"]][[name.expr]]$rp$bootstrap_results$annotation
    tmp.value <- list.gwas_inference[[name.gwas]][["inference"]][[name.expr]]$rp$bootstrap_results$bt_value # pvalue, 
    # rp$bootstrap_results : bt_value, pvalue
    # rp$full_results$parameters # gamma-estimates (is a named character vector)
    df.tmp <- data.frame(t(tmp.value)) # returns one-row data frame.
    colnames(df.tmp) <- tmp.annotations # set colnames
    ### Save in list
    list.dfs_meta_expr[[name.expr]][[name.gwas]] <- df.tmp
  }
}
### Combine
list.dfs_meta_expr <- lapply(list.dfs_meta_expr, bind_rows, .id="gwas")
names(list.dfs_meta_expr)



# ======================================================================================================= #
# ============================================= META PLOTS ============================================= #
# ======================================================================================================= #



for (name.expr in names(list.dfs_meta_expr)) {
  print(name.expr)
  
  # df <- list.dfs_meta_expr[["gtex.sub_tissue.avg_expr"]] # test
  df <- list.dfs_meta_expr[[name.expr]]
  
  df.gather <- df %>% gather(key="key", value="value", -gwas) # returns three column df (gwas, key, value)
  df.gather <- df.gather %>% group_by(gwas) %>% mutate(
    rank = base::rank(-value),
    p_value = stats::pnorm(value, lower.tail = F), # this is how the RolyPoly p-value is calculated
    significance = ifelse(p_value < 0.05, "*", ""), # asterisk or empty string
    rank_txt = paste0(rank, significance)
    ) 
  # ^ rank() ranks in ascending order (smallest number == higest rank). By taking the negative value, we get the descending order (largest value == highest rank)
  # ^ NB: large positive RolyPoly t-values correspond to significant p-values
  p <- ggplot(df.gather, aes(x=gwas, y=key, fill=value)) + geom_tile() + geom_text(aes(label=rank_txt))
  p <- p + labs(title=name.expr) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p <- p + scale_fill_gradient2()
  file.out <- sprintf("meta_params.heatmap.%s.pdf", name.expr)
  ggsave(p, filename=file.out, width=14, height=25)
  
  df.t <- t(df %>% column_to_rownames(var="gwas")) %>% # setting rownames and transposing. returns matrix with gwas as columnnames
    as.data.frame() # not really needed
  
  # p <- ggpairs(df.t) + labs(title=name.expr) #  # there is apparantly a small bug so ggpairs() only works if the columns have 'clean names'
  # file.out <- sprintf("meta_params.pairs.%s.pdf", name.expr)
  # ggsave(p, filename=file.out, width=30, height=30)
  
  p <- ggcorr(df.t, geom = "circle", hjust = 1, size = 2, layout.exp = 5) # label = TRUE # REF: ggcor: https://briatte.github.io/ggcorr/
  p <- p + labs(title=name.expr)
  file.out <- sprintf("meta_params.corr.%s.pdf", name.expr)
  ggsave(p, filename=file.out, width=12, height=12)
  
}


# ======================================================================================================= #
# ============================================= ROLYPOLY PLOTS ============================================= #
# ======================================================================================================= #

for (name.gwas in names(list.gwas_inference)) {
  print(sprintf("Processing GWAS=%s", name.gwas))
  
  dir.out <- sprintf("export.%s", name.gwas) # no trailing slash
  if (dir.exists(dir.out)) { # remove any existing folder
    print("Deleting previous output dir and creating a new")
    unlink(dir.out)
    dir.create(dir.out,showWarnings=F)
  } else {
    dir.create(dir.out,showWarnings=F)
  }
  
  for (name.expr in names(list.gwas_inference[[name.gwas]][["inference"]])) {
    print(sprintf("Processing EXPRESSION=%s", name.expr))
  
    ### Get rp object
    rp <- list.gwas_inference[[name.gwas]][["inference"]][[name.expr]]$rp
    # rp <- list.gwas_inference[["body_BMI_Locke2015.gene.10kb.squared.protein_coding_only"]][["inference"]][["gtex.tissue.gene_tpm.avg_expr"]]$rp
    
    ### tmp HACK (only needed when rp$raw_block_data has been removed)
    tmp.names <- names(rp$full_results$parameters)
    annotations <- tmp.names[!tmp.names %in% c("intercept")]
    rp$raw_block_data <- data.frame(matrix(ncol = length(annotations), nrow = 0))
    colnames(rp$raw_block_data) <- annotations

    ### Plot the log10(p) to rank tissues by the strength of their association
    p.pval_ranking <- plot_rolypoly_annotation_ranking(rp) + labs(title=name.expr) + theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=50.5, unit="points"))
    file.out.plot <- sprintf("%s/%s.plot_pval_ranking.%s.pdf", dir.out, name.gwas, name.expr)
    ggsave(p.pval_ranking, filename=file.out.plot, width=30 , height=10)
  
    ### Plot the TOP 10 log10(p) to rank tissues by the strength of their association
    top10_annotations <- rp$bootstrap_results %>%
      filter(!annotation %in% c("intercept", "maf_poly_1")) %>%
      arrange(-bt_value) %>%
      slice(1:100) %>%
      pull(annotation) %>%
      as.character()
    p.pval_ranking <- plot_rolypoly_annotation_ranking(rp) + labs(title=name.expr) + theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=50.5, unit="points"))
    p.pval_ranking <- p.pval_ranking + xlim(top10_annotations) + ylim(NA, -log10(rp$bootstrap_results$bp_value[rp$bootstrap_results$annotation==top10_annotations[1]]))
    p.pval_ranking
    file.out.plot <- sprintf("%s/%s.plot_pval_top10_ranking.%s.pdf", dir.out, name.gwas, name.expr)
    ggsave(p.pval_ranking, filename=file.out.plot, width=14, height=7)
  
    ### Plot gamma estimate with 95% confidence intervals
    p.gamma_est <- plot_rolypoly_annotation_estimates(rp) + labs(title=name.expr, x="gamma_hat_k") + theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=50.5, unit="points"))
    file.out.plot <- sprintf("%s/%s.plot_gamma_estimates.%s.pdf", dir.out, name.gwas, name.expr)
    ggsave(p.gamma_est, filename=file.out.plot, width=30, height=10)
  
    ### Export P-values
    file.out.table <- sprintf("%s/%s.table_pval_ranking.%s.csv", dir.out, name.gwas, name.expr)
    write_csv(rp$bootstrap_results %>% arrange(-bt_value), path=file.out.table)
  
    ### Export gamma-estimates
    df.gamma <- rp$full_results$parameters %>% # "rp$full_results$parameters" is a named character vector
      as_tibble() %>% # as_tibble/as.tibble/as_data_frame are aliases.
      rownames_to_column(var = "annotation") %>% # for some reason, you must call rownames_to_column() before calling arrange() - else you loose the rownames.
      arrange(desc(value))
    file.out.table <- sprintf("%s/%s.table_gamma.%s.csv", dir.out, name.gwas, name.expr)
    write_csv(df.gamma, path=file.out.table)
  }
}

# ======================================================================================================= #
# -============================================== PLOTS =============================================== #
# ======================================================================================================= #



# ======================================================================================================= #
# -============================================== XXXXXXX =============================================== #
# ======================================================================================================= #



# ======================================================================================================= #
# -============================================== XXXXXXX =============================================== #
# ======================================================================================================= #



# ======================================================================================================= #
# -============================================== LEFTOVERS =============================================== #
# ======================================================================================================= #

### For debugging
# name.gwas <- "body_BMI_Locke2015.gene.10kb.squared.protein_coding_only"
# name.expr <- "gtex.tissue.gene_tpm.avg_expr"
# tmp.annotations <- list.gwas_inference[[name.gwas]][["inference"]][[name.expr]]$rp$bootstrap_results$annotation
# tmp.value <- list.gwas_inference[[name.gwas]][["inference"]][[name.expr]]$rp$bootstrap_results$bt_value
# df.tmp <- data.frame(t(tmp.value)) # returns one-row data frame.
# colnames(df.tmp) <- tmp.annotations
# df.tmp



#### Extracting meta data
# magrittr, set_names
# 
# list.dfs_meta_expr <- list()
# for (name.gwas in names(list.gwas_inference)) {
#   for (name.expr in names(list.gwas_inference[[name.gwas]][["inference"]])) {
#     list.dfs_meta_expr[[name.expr]] <- set_names(list(list.gwas_inference[[name.gwas]][["inference"]][[name.expr]]$rp$bootstrap_results$bt_value), name.gwas) # using magrittr set_names()
#     # list.dfs_meta_expr[[name.expr]] <- list( UQ(rlang::sym(name.gwas)) = list.gwas_inference[[name.gwas]][["inference"]][[name.expr]]$rp$bootstrap_results$bt_value) # does not work: UQ(rlang::sym(name.gwas))
#   }
# }
# list.dfs_meta_expr
# 
# df.x <- lapply(list.dfs_meta_expr, function(l) {lapply(l, bind_rows, .id="gwas")})
