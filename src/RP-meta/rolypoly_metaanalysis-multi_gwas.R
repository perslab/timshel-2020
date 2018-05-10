############### SYNOPSIS ###################
# Analyzing RolyPoly results across multiple GWAS.

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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================================================= #
# =========================================== LOAD ROLYPOLY ============================================ #
# ======================================================================================================= #


### Get list of files
dir.rolypoly_objs <- "/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/out.rolypoly_objs-v1.clean/"
list.rolypoly_obj_files <- file.path(dir.rolypoly_objs, list.files(dir.rolypoly_objs, pattern="*.inference.RData")) # absolute path
head(list.rolypoly_obj_files)
print(sprintf("Got %s files to load", length(list.rolypoly_obj_files)))

### LOOP
list.gwas_inference <- list() # initialyzing empty list
for (file.rp in list.rolypoly_obj_files) {
  name.gwas <- stringr::str_split(basename(file.rp), "\\.")[[1]][3] # "out.rolypoly_objs.body_BMI_Locke2015.squared_tss_10kb.final.RData" --->  name.gwas = "body_BMI_Locke2015"
  print(sprintf("loading GWAS = %s ...", name.gwas))
  tmp.loadtime <- system.time(load(file.rp)) # list.rp_inference, list.run_parameters
  print(tmp.loadtime)
  
  list.gwas_inference[[name.gwas]] <- list("inference"=list.rp_inference, "parameters"=list.run_parameters, "gwas"=name.gwas)
}



# list.gwas_inference[[<GWAS>]][["inference"]][[<EXPR_DATASET>]][["rp"]]


# ======================================================================================================= #
# ====================================== SET NAMES (expression data) ==================================== #
# ======================================================================================================= #

names.list.rp_inference.orig <- names(list.rp_inference) # save orig names
names(list.rp_inference) <- str_replace_all(names.list.rp_inference.orig, c("^.*/(.*?)\\." = "\\1.", 
                                                                            "celltype_expr." = "", 
                                                                            ".hsapiens_orthologs.csv.gz" = ""))
print(names(list.rp_inference))

# ======================================================================================================= #
# -============================================== SET RP =============================================== #
# ======================================================================================================= #


### Things to extract
# rp$block_heritability_contribution (h_k^annot?)


for (name.expr in names(list.rp_inference)) {
  dir.out <- sprintf("export.%s", name.gwas) # no trailing slash
  if (dir.exists(dir.out)) { # remove any existing folder
    print("Deleting previous output dir and creating a new")
    unlink(dir.out)
    dir.create(dir.out,showWarnings=F)
  } else {
    dir.create(dir.out,showWarnings=F)
  }

  print(sprintf("Processing %s", name.expr))

  ### Get rp object
  rp <- list.rp_inference[[name.expr]]$rp
  
  ### tmp HACK (only needed when rp$raw_block_data has been removed)
  # tmp.names <- names(rp$full_results$parameters)
  # annotations <- tmp.names[!tmp.names %in% c("intercept")] 
  # rp$raw_block_data <- data.frame(matrix(ncol = length(annotations), nrow = 0))
  # colnames(rp$raw_block_data) <- annotations
  # rp$raw_block_data
  # 
  ### Plot the log10(p) to rank tissues by the strength of their association
  p.pval_ranking <- plot_rolypoly_annotation_ranking(rp) + labs(title=name.expr) + theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=50.5, unit="points"))
  file.out.plot <- sprintf("%s/%s.plot_pval_ranking.%s.pdf", dir.out, name.gwas, name.expr)
  ggsave(p.pval_ranking, filename=file.out.plot, width=30 , height=10)

  ### Plot the TOP 10 log10(p) to rank tissues by the strength of their association
  top10_annotations <- rp$bootstrap_results %>%
    filter(annotation != "intercept") %>%
    arrange(-bt_value) %>%
    slice(1:10) %>%
    pull(annotation) %>%
    as.character()
  p.pval_ranking <- plot_rolypoly_annotation_ranking(rp) + labs(title=name.expr) + theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=50.5, unit="points"))
  p.pval_ranking <- p.pval_ranking + xlim(top10_annotations)
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
# -============================================== XXXXXXX =============================================== #
# ======================================================================================================= #





