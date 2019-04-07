############### SYNOPSIS ###################
# Wrappers for running RolyPoly
# - 'robust' tryCatch() calls for RolyPoly GWAS-linking and making inference.
# - wrappers that loop over datasets

### TODO:
# load required libraries.
# pass all required arguments to functions in this script. Currently, functions in this script RELIES ON variables in the MAIN SCOPE to be available.

### OUTPUT: 
# ....

### REMARKS:
# ...

### REFERENCE:

############################################


# ======================================================================================================= #
# ========================================== Library ================================================== #
# ======================================================================================================= #

library(tidyverse)
library(doParallel)
library(foreach) # lib is also loaded as part of doParallel, but we just make sure to load it. Needed for "%dopar%"
library(rolypoly)

# ======================================================================================================= #
# ========================================== CONSTANTS ================================================== #
# ======================================================================================================= #

FLAG.SAVE_INFERENCE_PER_DATASET <- TRUE # save inference object after each dataset has been analyzed?
N_MAX_ANNOTATIONS <- Inf # Used for debugging/testing purpose. max number of annotations to analyze. Set to Inf to analyze all annotations in the dataset
# N_MAX_DATASETS <- Inf # Used for debugging/testing purpose. max number of data to analyze. Set to Inf to analyze all datasets

# ======================================================================================================= #
# ================================= UTILS (TODO: move to another file) ================================== #
# ======================================================================================================= #


get_runtime <- function(start_time) {
  # INPUT: start_time a "Sys.time()" object: class --> "POSIXct" "POSIXt"
  # REF: https://stackoverflow.com/questions/11017933/format-time-span-to-show-hours-minutes-seconds
  # OBS: there are also easier alternatives
  start_time <- as.POSIXct(start_time)
  # Sys.sleep(1) # for testing.
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  runtime_fmt_str <- format(.POSIXct(dt,tz="GMT"), "%H h:%M min:%S sec") # returns a character
  return(runtime_fmt_str)
}


# ======================================================================================================= #
# ================================= ROLYPOLY: PRE-LOAD GWAS FUNCTION ================================== #
# ======================================================================================================= #

run_rolypoly_precomputation.multi_attempt <- function(gwas_data,
                                                      snp_annotations,
                                                      block_annotation,
                                                      ld_folder,
                                                      n_attempts) {
  N_MAX_ATTEMPTS <- 10 # constant. max number of attempts
  n_attempts <- n_attempts + 1 # incrementing attempts
  
  rp <- tryCatch( # 'tryCatch()' will return the last evaluated expression  in case the "try" part was completed successfull
    {
      message(sprintf("Attempt = %s | *Try block* | Running GWAS precomputation attempt...", n_attempts))
      rolypoly <- list(); 
      class(rolypoly) <- 'rolypoly'
      rolypoly <- rolypoly:::rolypoly_load_gwas(rolypoly,
                                                gwas_data = gwas_data, 
                                                snp_annotations = snp_annotations, # character vector of column names in "gwas_data" data frame to include as SNP covariates in the model.
                                                # set snp_annotations to NULL if no covariates should be included.
                                                gwas_z_filter = -1, # gwas_z_filter is only 'active' if gwas_z_filter > 0
                                                add_spline = FALSE, # disable
                                                n_knots = 1, # only relevant if 'add_spline' == TRUE
                                                add_poly = TRUE, # add MAF 
                                                n_degree = 1)
      rolypoly <- rolypoly:::rolypoly_load_block_annotation(rolypoly, block_annotation) # load gene annotation
      rolypoly <- rolypoly:::rolypoly_link_blocks_and_gwas(rolypoly, ld_folder, run_parallel = T) # this takes quite a lot of time
      rolypoly$gwas <- NULL # now we can delete the GWAS slot.
      return(rolypoly) # You don't need to state the return value via return() but having 'rolypoly$gwas <- NULL' as last expression would make the return value "NULL".
    },
    error=function(cond) {
      message(sprintf("Attempt = %s | *Error block* | An error happend! Here's the original error message:", n_attempts))
      message(conditionMessage(cond)) # conditionMessage() needed for message to arrive when running in parallel mode. REF TryCatch messages for parallel R https://stackoverflow.com/a/45783896/6639640
      message() # gives extra space.
      if (n_attempts < N_MAX_ATTEMPTS) {
        message(sprintf("Attempt = %s | *Error block* | Will run function again...", n_attempts))
        ### Calling the function with the SAME arguments as it was called the first time.
        run_rolypoly_precomputation.multi_attempt(gwas_data,
                                                  snp_annotations,
                                                  block_annotation,
                                                  ld_folder,
                                                  n_attempts) # *OBS*: calling with updated variable
      } else {
        message(sprintf("Attempt = %s | *Error block* | REACHED MAX ATTEMPTS - will return 'NULL'", n_attempts))
        return(NULL) # return value if everything keeps failing
      }
    } # end error handling
  ) # end tryCatch
  return(rp) # return rolypoly
}


# ======================================================================================================= #
# =================================== ROLYPOLY INFERENCE FUNCTION ======================================= #
# ======================================================================================================= #

run_rolypoly_inference.multi_attempt <- function(rp.precomputed,
                                                 block_data,
                                                 do.parallel_bootstrap,
                                                 n_attempts) {
  N_MAX_ATTEMPTS <- 10 # constant. max number of attempts
  n_attempts <- n_attempts + 1 # incrementing attempts
  
  rp <- tryCatch(
    {
      # 'tryCatch()' will return the last evaluated expression  in case the "try" part was completed successfull
      message(sprintf("Attempt = %s | *Try block* | running RolyPoly inference...", n_attempts))
      
      tmp.time <- system.time(rolypoly <- rolypoly:::rolypoly_load_block_data(rolypoly = rp.precomputed, # RolyPoly object precomputed (linked GWAS and gene annotation data)
                                                                              block_data = block_data)) # load expression data
      print("time rolypoly_load_block_data:")
      print(tmp.time)
      tmp.time <- system.time(rolypoly <- rolypoly:::rolypoly_perform_inference(rolypoly = rolypoly, 
                                                                                bootstrap_iters = N_BOOTSTRAP_ITERS, # use at least 200
                                                                                outlier_threshold = -1, # -1 means inactive. "outlier_threshold threshold for performing robust regression, still experimental."
                                                                                run_light = T, # if we want to throw away bootstrap data, and save memory 
                                                                                run_parallel = do.parallel_bootstrap)) # run bootstrapping in parallel 
      print("time rolypoly_perform_inference:")
      print(tmp.time)
      rolypoly # return value. Important to have this as the last statement.
    },
    error=function(cond) {
      message(sprintf("Attempt = %s | *Error block* | An error happend! Here's the original error message:", n_attempts))
      # message(cond) # ---> error message does not arrive from parallel workers to 'master'
      message(conditionMessage(cond)) # conditionMessage() needed for message to arrive when running in parallel mode. REF TryCatch messages for parallel R https://stackoverflow.com/a/45783896/6639640
      message() # gives extra space.
      if (n_attempts < N_MAX_ATTEMPTS) {
        message(sprintf("Attempt = %s | *Error block* | Will run function again...", n_attempts))
        ### Calling the function with the SAME arguments as it was called the first time.
        run_rolypoly_inference.multi_attempt(rp.precomputed = rp.precomputed,
                                             block_data = block_data,
                                             do.parallel_bootstrap = do.parallel_bootstrap,
                                             n_attempts = n_attempts) # *OBS*: calling with updated variable
      } else {
        message(sprintf("Attempt = %s | *Error block* | REACHED MAX ATTEMPTS - will return 'NULL'", n_attempts))
        return(NULL) # return value if everything keeps failing
      }
    } # end error handling
  ) # end tryCatch
  return(rp) # return rolypoly
}




# ======================================================================================================= #
# ========================================== *PARALLEL ANNOTATION* ====================================== #
# ================================ Run RolyPoly over multiple expression data =========================== #
# ======================================================================================================= #


wrapper.run_rolypoly_univariate_over_expr_datasets <- function(list.df_expr, 
                                                               rp.precomputed, 
                                                               do.parallel_annotations=TRUE,
                                                               do.parallel_bootstrap=TRUE) {
  ### INPUT
  # list.df_expr              a named list of data frames with expression data sets. Names should be names of expression data.
  # rp.precomputed            gwas_linked rolypoly object.
  # do.parallel_annotations  boolean. If TRUE then run annotations in parallel.
  #                           we recommend setting to TRUE for univariate regressions, for the following reasons:
  #                             1. the bootstrapping runs really fast for small regressions and parallel bootstraping is much more overhead.
  #                             2. rolypoly_load_block_data() takes 1-2 min per annotation (depending on the number of genes).
  #                           [DISREGARD THIS --> ] if TRUE the bootstrapping *will not be run in parallel* because the 'parallel backend' is used on the parallel annotations.
  #                           [DISREGARD THIS --> ] if TRUE, the attempt to run bootstrap in parallel might cause a harmless warning: "Warning message: executing %dopar% sequentially: no parallel backend registered"
  # do.parallel_bootstrap   boolean. If TRUE then run bootstrap in parallel.
  #                         *OBS*: if do.parallel_bootstrap=T and do.parallel_annotations=T, then N_CORES x N_CORES processes will be created.
  
  
  if (is.null(names(list.df_expr))) { # we require names
    stop("Error: list.df_expr is not a *named* list. Fix this by adding names to the list.")
  }
  
  # determine parallelization mode
  if (do.parallel_annotations) {
    "%loop_function%" <- `%dopar%`
  } else {
    "%loop_function%" <- `%do%`
  }
  
  ### Start parallel
  ### We put the registerDoParallel within the for-loop, using that it will make the parallelization more stable
  print("Starting registerDoParallel...")
  # cl <- parallel::makeCluster(N_CORES, type="FORK")
  # registerDoParallel(cl) # The registerDoParallel function is used to register the parallel backend with the foreach package.
  registerDoParallel(N_CORES)
  print(sprintf("Number of registered workers (getDoParWorkers()): %s", getDoParWorkers()))
  ### R parallel computing and zombie processes
  # REF: https://stackoverflow.com/questions/25388139/r-parallel-computing-and-zombie-processes
  # REF: https://stackoverflow.com/questions/25348607/how-to-stop-r-from-leaving-zombie-processes-behind
  # 12367 12370 12367 35481 pts/34   12367 S+   132983994   0:00      \_ /bin/sh -c Rscript run_rolypoly_multi_expr_datasets-v3.R --gwas_name body_BMI_Yengo2018 --run_name tss.10kb.none.all_genes --outdir /projects/timshe
  # 12370 12374 12367 35481 pts/34   12367 R+   132983994 2194:08      |   \_ /usr/lib64/R/bin/exec/R --slave --no-restore --file=run_rolypoly_multi_expr_datasets-v3.R --args --gwas_name body_BMI_Yengo2018 --run_name tss.10kb.n
  # 12374 78262 12367 35481 pts/34   12367 Z+   132983994   3:35      |       \_ [R] <defunct>
  #   12374 78264 12367 35481 pts/34   12367 Z+   132983994   3:25      |       \_ [R] <defunct>

  list.rp_inference_full <- list() # initialyzing return list. 
  # ^ Will become a list of lists of lists: list.rp_inference_full[expr_data_name][name.annotation] --> list of rp, runtime
  
  n.expr_data_sets <- length(list.df_expr)
  counter_expr_data <- 0
  for (name.expr_data in names(list.df_expr)) { # looping over expression data sets
    start_time.expr_data <- Sys.time() # start timing
    
    counter_expr_data <- counter_expr_data + 1
    df.expr <- list.df_expr[[name.expr_data]] # get expression data frame
    print(sprintf("========== Dataset #%s/#%s | Status = START RUNNING | wrapper.run_rolypoly_over_expr_datasets | list element = %s ==========", counter_expr_data, n.expr_data_sets, name.expr_data)) # prints "character(0)" if list has no names (i.e. names(list) is NULL)
    # n_annotations <- ncol(df.expr)
    n_annotations <- min(ncol(df.expr), N_MAX_ANNOTATIONS) # for testing
    
    ## .packages = c("tidyverse", "rolypoly")
    ## TODO: set names of list --> .final = function(x) setNames(x, names(foo)) [+ .inorder=T] # REF: https://stackoverflow.com/questions/27276269/foreach-keep-names
    ## *OBS*: remember to use double squarebrackets for the list assignment or else the assigment will be *WRONG* (and you will get the "Warning message: number of items to replace is not a multiple of replacement length")
    list.rp_inference_full[[name.expr_data]] <- foreach (idx_col.annotation=1:n_annotations, # looping over annotations/cell_types
                                                    .inorder=FALSE # .inorder=FALSE can increase performance and we do not need the annotations to be in order.
                                                    ) %loop_function% { # ") %loop_function% {" should be on the same line, otherwise a syntax error is thrown
      name.annotation <- colnames(df.expr)[idx_col.annotation]
      print(sprintf("========== Annotation #%s/#%s | Status = Running inference | Annotation = %s ==========", idx_col.annotation, n_annotations, name.annotation))
      df.annotation <- df.expr %>% select(name.annotation) %>% as.data.frame()
      # print(head(df.annotation))
      runtime <- system.time( # system.time(): returns a names numerical vector of class "proc_time"
        rp <- run_rolypoly_inference.multi_attempt(rp.precomputed = rp.precomputed,
                                                   block_data = df.annotation, 
                                                   do.parallel_bootstrap = do.parallel_bootstrap, # *VARIABLE*
                                                   n_attempts = 0)# *OBS*: important to make first call with n_attempts=0
      )
      ###
      print(sprintf("annotation = %s | RUNTIME run_rolypoly_inference.multi_attempt:", name.annotation))
      print(runtime)
      
      ### Deleting data slots from RolyPoly object
      # We do this to reduce later load times and save space (storage and memory).
      rp$data <- "REMOVED" # *IMPORTANT* gene level information (GWAS-linked) [takes up 98% of the size of the object]
      rp$blocks <- "REMOVED" # gene==block annotation
      rp$raw_block_data <- "REMOVED" # gene==block expression data
      rp$full_results$model <- "REMOVED" # *IMPORTANT* this is like the 'keep_mode'=F in the main_wrapper.R [this slot use around ~400 MB]
      # rp$snp_annotations <- "REMOVED" # rp$snp_annotations is only a character vector?

      
      ### "return value"
      # list.rp_inference_full[[name.expr_data]][[name.annotation]] <- list("rp"=rp, "runtime"=runtime, "name.expr_data"=name.expr_data, "annotation"=name.annotation) # *OBS*: notice the list[[x]] indexing: list[x] is not enough
      list("rp"=rp, "runtime"=runtime, "name.expr_data"=name.expr_data, "name.annotation"=name.annotation) # *OBS*: notice the list[[x]] indexing: list[x] is not enough
      # *OBS* the last expression in this block is the 'return value', so DO NOT put any expressions after this.
    } # end foreach annotation
    print(sprintf("========== Dataset #%s/#%s | Status = DONE | wrapper.run_rolypoly_over_expr_datasets | list element = %s ==========", counter_expr_data, n.expr_data_sets, name.expr_data)) # prints "character(0)" if list has no names (i.e. names(list) is NULL)
    print(sprintf("Runtime for expression dataset: %s", get_runtime(start_time.expr_data)))
    print(sprintf("Size list.rp_inference_full: %s", format(object.size(list.rp_inference_full),units="MB")))
    
    
    ### save object for each dataset's finished inference
    if (FLAG.SAVE_INFERENCE_PER_DATASET) {
      print(names(list.rp_inference_full))
      print("temporary saving inference object...")
      list.rp_inference <- list()
      for (name.expr_data in names(list.rp_inference_full)) {
        print(name.expr_data)
        list.rp_inference[[name.expr_data]] <- list() # TODO: not sure this is needed. Try to remove it.
        for (i.anno in seq_along(list.rp_inference_full[[name.expr_data]])) { # *OBS*: looping over i.anno
          name.annotation <- list.rp_inference_full[[name.expr_data]][[i.anno]][["name.annotation"]]
          print(sprintf("%s=%s", i.anno, name.annotation))
          # file.out.inference.csv_bootstrap <- sprintf("%s.bootstrap.%s.%s.csv", file.out.inference.csv_basename, name.expr_data, name.annotation)
          # print("writing csv file 1...")
          # write_csv(df.bootstrap_results, path = file.out.inference.csv_bootstrap)
          # file.out.inference.csv_gamma <- sprintf("%s.gamma.%s.%s.csv", file.out.inference.csv_basename, name.expr_data, name.annotation)
          # print("writing csv file 2...")
          #write_csv(df.gamma, path = file.out.inference.csv_gamma)
          # h_j_gene == "rolypoly$block_values"
          # h_k_annot == rolypoly$block_heritability_contribution
          list.rp_inference[[name.expr_data]][[name.annotation]] <- list("bootstrap_results"=list.rp_inference_full[[name.expr_data]][[i.anno]][["rp"]]$bootstrap_results, # p-values enirchment [data.frame: annotations x 'estimates']. Size per anno=insignificant
                                                                "gamma"=list.rp_inference_full[[name.expr_data]][[i.anno]][["rp"]]$full_results$parameters, # gamma estimates [named numerical: n=annotations]. Size per anno=insignificant
                                                                "block_values"=list.rp_inference_full[[name.expr_data]][[i.anno]][["rp"]]$block_values, # h_j_gene [named numerical: n=genes]. Size per anno=[GTEx_49126_genes="3.7 Mb"; kme_5000_genes="0.3 Mb"]
                                                                "block_heritability_contribution"=list.rp_inference_full[[name.expr_data]][[i.anno]][["rp"]]$block_heritability_contribution, # h_k_annot [scalar]. Size per anno=insignificant
                                                                "expected_block_values"=list.rp_inference_full[[name.expr_data]][[i.anno]][["rp"]]$expected_block_values) # XXX [data.table: genes x 3 columns (g_hat, g, label)]. Size per anno=[GTEx_49126_genes="3.7 Mb"; kme_5000_genes="0.3 Mb"]
                                                                
        } #  save, end foreach anno
      } # save, end foreach dataset
      
      ### Writing list.rp_inference
      runtime <- system.time(save(list.rp_inference, list.run_parameters, file=file.out.inference.tmp_red))
      print("Time saving inference object [list.rp_inference]:")
      print(runtime)
      
      ### Saving temporary list.rp_inference_full containing 'full' RolyPoly objects --> WRITING FILE IS UNBELIEVABLE SLOW!! Will practically never finish.
      # runtime <- system.time(save(list.rp_inference_full, list.run_parameters, file=file.out.inference.tmp))
      # print("Time saving inference object [list.rp_inference_full]:")
      # print(runtime)
      
    } # end if FLAG.SAVE_INFERENCE_PER_DATASET
    
    # print("Calling garbage collection...")
    # gc()
    
  } # end foreach data set

  
  ### stop cluster
  print("Stopping cluster...")
  # stopCluster(cl)
  stopImplicitCluster() # stop cluster [stopCluster() should not be used or?]
  print(sprintf("After stopping cluster: number of registered workers (getDoParWorkers()): %s", getDoParWorkers()))
  
  # return(list.rp_inference_full)
  return(list.rp_inference) # return 'reduced' list (selectively extracted RolyPoly elements)
}


### DEBUGGING rolypoly_load_block_data [keep for now]
# list.rp_inference.run <- lapply(seq_along(list.df_expr.run), wrapper.run_rolypoly_univariate_over_expr_datasets, list=list.df_expr.run, rp.precomputed=rp.gwas_linked) # "list" is here the extra (named) argument to wrapper.run_rolypoly_over_expr_datasets()
# trace("rolypoly_load_block_data", browser, where=asNamespace("rolypoly"))
# trace("rolypoly_load_block_data", edit=T, where=asNamespace("rolypoly")) # where is needed because 'rolypoly_load_block_data' is not part of the package 'exported' functions
# untrace("rolypoly_load_block_data", where=asNamespace("rolypoly"))
# rolypoly:::rolypoly_perform_inference


