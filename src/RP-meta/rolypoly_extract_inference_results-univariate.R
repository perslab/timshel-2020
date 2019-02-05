############### SYNOPSIS ###################
# Analyzing multiple RolyPoly run

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################


# rstudioapi::getActiveDocumentContext()


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

# suppressMessages(library(rolypoly))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))

suppressMessages(library(doParallel))
suppressMessages(library(foreach))

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================================================= #
# ============================================== GET FILES ============================================== #
# ======================================================================================================= #

### LDSC munge GWAS | new 2019-01-30
DIR.DATA_INFERENCE <- "/scratch/rolypoly_out_scratch/out.rolypoly_objs-v3.univariate/"
DIR_OUT.PREFIX <- "inference_rp_ldscmunge"
DIR_OUT.COMBINED_TAR <- "export-combined.rp.v3"
list.inference_files <- list.files(path=DIR.DATA_INFERENCE,  pattern="*.inference.RData")

### Main
# DIR.DATA_INFERENCE <- "/scratch/sc-genetics/out.rolypoly_objs-v3.univariate-nboot100/"
# DIR_OUT.PREFIX <- "inference_rp"
# DIR_OUT.COMBINED_TAR <- "export-combined.v3.nboot100"
# list.inference_files <- list.files(path=DIR.DATA_INFERENCE,  pattern="*.inference.RData")

### Main (HM3) - used before 190130
# DIR.DATA_INFERENCE <- "/scratch/sc-genetics/out.rolypoly_objs-v3.hm3_protein_coding_only.univariate-nboot100/"
# DIR_OUT.PREFIX <- "inference_rp_hm3"
# DIR_OUT.COMBINED_TAR <- "export-combined.rp_hm3.v3.nboot100"
# list.inference_files <- list.files(path=DIR.DATA_INFERENCE,  pattern="*.inference.RData")

### LDSC gwas
# DIR.DATA_INFERENCE <- "/scratch/sc-genetics/out.rolypoly_objs-v3.ldsc.univariate-nboot100/"
# DIR_OUT.PREFIX <- "inference_ldsc"
# DIR_OUT.COMBINED_TAR <- "export-combined.ldsc.v3.nboot100"
# list.inference_files <- list.files(path=DIR.DATA_INFERENCE,  pattern="*.inference.RData")

### NULL gwas
# DIR.DATA_INFERENCE <- "/scratch/sc-genetics/out.rolypoly_objs-v3.NULL_GWAS/"
# DIR_OUT.PREFIX <- "inference_null_gwas"
# DIR_OUT.COMBINED_TAR <- "export-combined.null_gwas.v3.nboot100"
# list.inference_files <- list.files(path=DIR.DATA_INFERENCE,  pattern="*.inference.tmp_red.RData") # NULL GWAS




### Source annotations
source("/raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/constants-annotation_name_mapping.R")

do.parallel_annotations <- TRUE
# do.parallel_annotations <- FALSE
N_CORES <- 50



# ======================================================================================================= #
#============================================= FUNCTION ================================================= #
# ======================================================================================================= #

extract_data <- function(list.rp_inference) {
  
  # list.rp_inference <- list.rp_inference # TESTING
  # name.expr_data <- "gtex.sub_tissue_gene_tpm.avg_expr" # TESTING
  # name.anno <- "Adipose - Subcutaneous" # TESTING
  # obj_type="full_rp"
  
  list.bootstrap <- list()
  list.gamma <- list()
  for (name.expr_data in names(list.rp_inference)) {
    print(name.expr_data)
    list.bootstrap[[name.expr_data]] <- list()
    list.gamma[[name.expr_data]] <- list()
    for (name.anno in names(list.rp_inference[[name.expr_data]])) {
      print(sprintf("annot=%s", name.anno))
      list.bootstrap[[name.expr_data]][[name.anno]] <- list.rp_inference[[name.expr_data]][[name.anno]][["bootstrap_results"]][3,] # 1 row data.frame with the annotation results
      list.gamma[[name.expr_data]][[name.anno]] <- list.rp_inference[[name.expr_data]][[name.anno]][["gamma"]] %>% # named numeric
        as_tibble() %>% # as_tibble/as.tibble/as_data_frame are aliases.
        rownames_to_column(var = "annotation") %>% # for some reason, you must call rownames_to_column() before calling arrange() - else you loose the rownames.
        spread(key=annotation, value=value) %>% # returns a single-row dataframe
        select(intercept, maf_poly_1, everything()) # making the 'annotation' column go last
      colnames(list.gamma[[name.expr_data]][[name.anno]]) <- c("intercept", "maf_poly_1", "gamma")
      # rename(annotation := UQ(rlang::sym(name.anno)) ) # does not work because e.g. "Adipose - Subcutaneous" is transformed into Adipose...Subcutaneous
    } # end foreach anno
    list.bootstrap[[name.expr_data]] <- bind_rows(list.bootstrap[[name.expr_data]]) %>% arrange(annotation)
    list.gamma[[name.expr_data]] <- bind_rows(list.gamma[[name.expr_data]], .id="annotation") %>% arrange(annotation)
  }
  list.return <- list("bootstrap"=list.bootstrap, "gamma"=list.gamma)
  return(list.return)
}


add_annotations <- function(df, name.expr_data) {
  ### OUTPUT: adds columns 'category' and 'name_clean' to df
  if (grepl("gtex", name.expr_data)) {
    df.category <- df.category.gtex
  } else if (grepl("maca", name.expr_data)) {
    df.category <- df.category.maca
  } else if (grepl("depict", name.expr_data)) {
    df.category <- df.category.depict
  } else if (grepl("mousebrain", name.expr_data)) {
    df.category <- df.category.mousebrain
  # } else if (grepl("tabula_muris", name.expr_data)) {
  #   df.category <- df.category.tabula_muris
  #  ---> NOT DONE YET
  } else {
    print("WARNING: add_annotations() got unknown pattern for name.expr_data. Will do dummy operation. Plot might not have meaningful coloring.")
    df.category <- df %>% transmute(name_r=annotation, name_clean=annotation) # 'dummy' operation, but to make it work for ANY dataset (i.e that does not match the 'grepl' pattern)
    # stop("Got unexpected name for name.expr_data")
  }
  df.res <- left_join(df, df.category, by=c("annotation"="name_r"))
  return(df.res)
}

  

# ======================================================================================================= #
# =============================================== MAIN ================================================= #
# ======================================================================================================= #


### Determine parallelization mode
# REF - read this about GUI/Rstudio and parallel: https://stackoverflow.com/a/25126473/6639640
if (do.parallel_annotations) {
  print("Starting registerDoParallel...")
  # cl <- parallel::makeCluster(N_CORES, type="FORK")
  # registerDoParallel(cl) # The registerDoParallel function is used to register the parallel backend with the foreach package.
  registerDoParallel(N_CORES)
  print(sprintf("Number of registered workers (getDoParWorkers()): %s", getDoParWorkers()))
  "%loop_function%" <- `%dopar%`
} else {
  "%loop_function%" <- `%do%`
}

list.dummy <- foreach (i=seq_along(list.inference_files)) %loop_function% {
  ### DEBEGGUNG
  # i <- 1
  
  file.rp <- list.inference_files[[i]]
  print(file.rp)
  load(sprintf("%s/%s", DIR.DATA_INFERENCE, file.rp))
  list.res <- extract_data(list.rp_inference) # **CALL FUNCTION TO EXTRACT DATA**
  list.bootstrap <- list.res[["bootstrap"]]
  
  name.gwas_run <- stringr::str_replace_all(basename(file.rp), c("rolypoly_objs."="", ".inference.RData"=""))
  print(name.gwas_run)
  
  ### GWAS_main_out_dir
  dir.out <- sprintf("%s.%s", DIR_OUT.PREFIX, name.gwas_run) # no trailing slash
  if (dir.exists(dir.out)) { # remove any existing folder
    print("Deleting previous output dir and creating a new")
    unlink(dir.out)
    dir.create(dir.out,showWarnings=F)
  } else {
    dir.create(dir.out,showWarnings=F)
  }
  
  
  ### BOOTSTRAP
  list.bootstrap.all_expr_data <- list()
  for (name.expr_data in names(list.bootstrap)) {
    print(name.expr_data)
    df <- list.bootstrap[[name.expr_data]] # data frame with bootstrap results

    
    ### EXPR_main_out_dir (TMP)
    dir.out.expr_collection <- sprintf("%s.0.%s", DIR_OUT.PREFIX, name.expr_data) # no trailing slash
    if (!dir.exists(dir.out.expr_collection)) { # remove any existing folder
      dir.create(dir.out.expr_collection,showWarnings=F)
    }
    
    
    ### MODIFY DF
    df <- add_annotations(df, name.expr_data) # adds columns 'category' and 'name_clean'
    fdr_cutoff <- 0.05/nrow(df)
    df$fdr_significant <- if_else(df$bp_value <= fdr_cutoff, TRUE, FALSE)

    ### saving df
    list.bootstrap.all_expr_data[[name.expr_data]] <- df # saving
    
    
    file.out.csv <- sprintf("%s/table.pvals.%s.csv", dir.out, name.expr_data)
    write_csv(df %>% arrange(bp_value), path=file.out.csv)
    write_csv(df %>% arrange(bp_value), path=sprintf("%s/table.pvals.%s.csv", dir.out.expr_collection, name.gwas_run)) # *** TMP***
    
    ### Prepating for barplot: ORDERING name_clean by the 'category'
    df <- df %>% arrange(category) %>% mutate(name_clean = factor(name_clean, unique(name_clean)))
    # The key to (re)ordering factor levels, is that it MUST be a unique set of values supplied to 'levels'
    # ALTERNATIVE #2 - I COULD NOT FIND A WAY TO MAKE forcats::fct_reorder work by ordering a character column
    # fct_reorder(f, factor(df.lab.long$group), fun=XXX)
    
    
    map_size <- c("TRUE"=6, "FALSE"=3)
    
    #ggplot(df, aes(x=annotation, y=-log10(bp_value))) + 
    ggplot(df, aes(x=name_clean, y=-log10(bp_value))) + 
      geom_point(aes(color=category, size=fdr_significant)) +
      geom_hline(yintercept = -log10(fdr_cutoff), color="red") +
      theme(axis.text.x =  element_text(angle = 45, hjust = 1)) +
      theme(plot.margin=margin(t=5.5, r=5.5, b=5.5, l=60.5, unit="points")) +
      labs(title=name.expr_data, x="Annotation", y="-log10(P-val)") + 
      scale_size_manual(name="DUMMY", values=map_size) + 
      guides(size=FALSE) # hide size legend
    file.out.plot <- sprintf("%s/plot.pvals.%s.pdf", dir.out, name.expr_data)
    ggsave(filename = file.out.plot, width= 35, height = 10)
    ggsave(filename = sprintf("%s/plot.pvals.%s.pdf", dir.out.expr_collection, name.gwas_run), width= 35, height = 10) # *** TMP***
    # break
  }
  ### combine all results
  df.bootstrap.all_expr_data <- bind_rows(list.bootstrap.all_expr_data, .id="dataset")
  file.out.csv <- sprintf("%s/table.pvals.%s.csv", dir.out, "ALL_DATA")
  write_csv(df.bootstrap.all_expr_data %>% arrange(bp_value), path=file.out.csv)
  
  # ### make hypothalamus plot
  # df.bootstrap.all_expr_data %>% 
  #   filter(dataset %in% c("hypothalamus.campbell.avg_expr", 
  #                         "hypothalamus.chen.avg_expr", 
  #                         "hypothalamus.lira.avg_expr", 
  #                         "hypothalamus.romanov.avg_expr")) %>%
  #   ggplot(aes(x=annotation, y=-log10(bp_value))) + 
  #   geom_col() +
  #   geom_hline(yintercept = -log10(0.05), color="red") +
  #   theme(axis.text.x =  element_text(angle = 45, hjust = 1)) +
  #   labs(title="All hypothalamus data (analyzed individually)")
  # file.out.plot <- sprintf("%s/plot.pvals.%s.pdf", dir.out, "ALL_HYPO_ANALYZED_INDIVIDUALLY")
  # ggsave(filename = file.out.plot, width= 25, height = 10)
  
}
  

list.dummy[[1]]


##### EXTRACT DATA #######
cmd0 <- sprintf("rm -r %s %s.tar.gz", DIR_OUT.COMBINED_TAR, DIR_OUT.COMBINED_TAR) # REMOVE ANY OLD FOLDER. This is ok - we can always regerate it
system(cmd0)

cmd <- sprintf("mkdir -p %s", DIR_OUT.COMBINED_TAR) # make dir
system(cmd)
cmd1 <- sprintf("mv %s* %s/", DIR_OUT.PREFIX, DIR_OUT.COMBINED_TAR) # move dirs
system(cmd1)
cmd2 <- sprintf("tar -czvf %s.tar.gz %s", DIR_OUT.COMBINED_TAR, DIR_OUT.COMBINED_TAR)
system(cmd2)
dir()

print(sprintf("time scp ygg:%s/%s.tar.gz ~/Downloads/", getwd(), DIR_OUT.COMBINED_TAR))


# DIR.export-combined.v3.nboot100
# mkdir $DIR_OUT; mv export.* $DIR_OUT/; tar -czvf $DIR_OUT.tar.gz $DIR_OUT; ll












