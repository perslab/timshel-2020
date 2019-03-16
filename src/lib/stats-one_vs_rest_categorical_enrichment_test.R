
library(tidyverse)

library(doParallel)
library(foreach) # lib is also loaded as part of doParallel, but we just make sure to load it. Needed for "%dopar%"


# ======================================================================= #
# ======================== ENRICHMENT TESTS ============================= #
# ======================================================================= #

do_enrichment_test_category <- function(df, category.select, colname_test_for_enrichment, colname_statistic, do_empirical_pvalue) {
  ### DOCS: function to perform enrichment test
  ### The function implements a "category.select vs rest" one-sided wilcox.test [alternative="greater"].
  ### POTENTIAL LIMITATION OF TEST: 
  #   Cell-types within a category are not independent of each other. They are likely to be correlated. 
  #   This may violate the wilcox test and is not easily solved by the empirical p-value.
  
  ### INPUT
  # df:                           data frame.
  # colname_test_for_enrichment:  character. columnname in df that defines the categories.
  # category.select:              character. 'selected' category to compute enrichment statistic for.
  # colname_statistic:            character. columname in df that defines the 'value' that the categories will be compared based on.
  ### OUTPUT
  # pvalue:                       a single numeric value.
  
  print(category.select)
  df.now <- df %>% mutate(category.tmp = if_else(!!rlang::sym(colname_test_for_enrichment)==category.select, category.select, "Other"))
  x <- df.now %>% filter(category.tmp==category.select) %>% pull(!!rlang::sym(colname_statistic))
  y <- df.now %>% filter(category.tmp=="Other") %>% pull(!!rlang::sym(colname_statistic))
  wtest_obs <- wilcox.test(x,y , alternative="greater") # alternative hypothesis: x > y.
  # wtest_obs <- wilcox.test(bt_value ~ category.tmp, data = df.now, alternative="greater") # *OBS*: using formula is DANGEROUS (when using a one-sided test) because the levels in 'category.tmp' are ordered alphabetically
  
  ### Empirical p-value by randomizing category.tmp N_NULL times.
  if (do_empirical_pvalue) {
    set.seed(0)
    N_NULL <- 1000
    vector_null <- vector(mode="double", length=N_NULL)
    for (i in seq(N_NULL)) {
      print(sprintf("Null %s", i))
      df.now$category.tmp <- sample(df.now$category.tmp, replace=FALSE) # randomize assignments (permute)
      # *IMPORTANT*: to account for the 'category size', we keep the same number of "category.select" and "Other" as in the 'oberserved data'
      # *OBS*: A limitation of this sampling is that we do not keep the correlation structure in our data.
      #      ^ We know that cell-types within a category are likely to be correlated (e.g. cell-types from the same region) and hence their value in colname_statistic will be correlated.
      #      ^ This effectively means that cell-types and their colname_statistic might be considered 'blocks' of values that go together.
      # ^ https://stackoverflow.com/questions/13765972/how-to-randomize-a-vector
      x <- df.now %>% filter(category.tmp==category.select) %>% pull(!!rlang::sym(colname_statistic))
      y <- df.now %>% filter(category.tmp=="Other") %>% pull(!!rlang::sym(colname_statistic))
      vector_null[i] <- wilcox.test(x,y , alternative="greater")$statistic
    }
    idx.insertion <- findInterval(wtest_obs$statistic, sort(vector_null)) # find the rank of obs inside the sorted null.
    # return values of findInterval(x, vec) can take on values of {0...length(vec)}
    n_null_gt_obs <- length(vector_null)-idx.insertion # number of null observations *greather than* our observed value. 
    # ^ e.g. if an element in idx.insertion was the most extreme in the null, we have length(vector_null) == idx.insertion so n_null_gt_obs=0 for that element
    pval_empirical <- (n_null_gt_obs+1)/(length(vector_null)+1)
    # *OBS*: Note that the above formula includes a commonly used pseudocount (North, Curtis, and Sham, 2002; Knijnenburg et al., 2009) to avoid P-values of zero.
    ### ALTERNATIVE shorter but less transparant version: 
    # pval_empirical <- 1-(idx.insertion/(length(vector_null)+1)) # *OBS*: here you should not have +1 to idx.insertion because idx.insertion==length(vector_null) if our observed value is the most extreme
  }
  if (do_empirical_pvalue) {
    return(pval_empirical) 
  } else {
    return(wtest_obs$p.value)
  }
}


run_enrichment <- function(df.dataset, colname_test_for_enrichment, colname_statistic, do_empirical_pvalue, run_parallel=F, n_cores=20) {
  ### DOCS: WRAPPER function to run enrichment tests
  
  ### Determine parallelization mode
  # REF - read this about GUI/Rstudio and parallel: https://stackoverflow.com/a/25126473/6639640
  if (run_parallel) {
    print("Running in parallel mode")
    registerDoParallel(n_cores)
    print(sprintf("Number of registered workers (getDoParWorkers()): %s", getDoParWorkers()))
    "%loop_function%" <- `%dopar%`
  } else {
    "%loop_function%" <- `%do%`
  }
  
  ### iter_control is the unique values (categories) in colname_test_for_enrichment
  iter_control <- df.dataset %>% 
    distinct(!!sym(colname_test_for_enrichment)) %>%
    pull()
  # iter_control <- iter_control[1:10] # DEV
  list.res <- foreach (idx.foreach=1:length(iter_control), # columns
                       # .export(df.dataset, colname_test_for_enrichment, colname_statistic),
                       .inorder=TRUE, # .inorder=FALSE can increase performance, but we want things to be in order so we can set the names
                       .packages=c("tidyverse")) %dopar% {
                         name_element.foreach <- iter_control[idx.foreach] # an element in iter_control
                         print(paste("processing #", idx.foreach, "/#", length(iter_control), sep="" ))
                         print(paste("name_element.foreach is:", name_element.foreach))
                         
                         pvalue <- do_enrichment_test_category(df=df.dataset, 
                                                               category.select=name_element.foreach, # SELECTOR
                                                               colname_test_for_enrichment=colname_test_for_enrichment, 
                                                               colname_statistic=colname_statistic,
                                                               do_empirical_pvalue=do_empirical_pvalue)
                         df.return <- tibble(enrichment_val=pvalue) # we return a tibble to enable combining results via bind_rows()
                       }
  
  
  ### Setting names of result list
  # Since the ".inorder" is TRUE the list is combined in the same order that the elements were submitted.
  names(list.res) <- iter_control
  ### Combine results into tibble
  df.wtest_res <- bind_rows(list.res, .id="category")
  # ^ gives a two-column tibble: category and enrichment_val.
  
  ### Count number of observations in each group
  df.n <- df.dataset %>% 
    count(!!sym(colname_test_for_enrichment)) %>%
    rename(n_obs = n)
  ### Combine and order
  df.wtest_res <- df.wtest_res %>% 
    left_join(df.n, by=c(category=colname_test_for_enrichment)) %>%
    arrange(enrichment_val)
  
  return(df.wtest_res)
  ### Example return value
  #   category               enrichment_val n_obs
  # 1 Cortex                       0.000999    33
  # 2 Hippocampus                  0.000999    21
}



# ======================================================================= #
# =================== OLD CODE USING sapply as wrapper ================== #
# ======================================================================= #
#### OK TO DELETE!


# 
# ### ===== Select 'dataset' =====
# # Generic (neurons only)
# df.dataset <- df.ldsc_cts %>% filter(Class=="Neurons")
# colname_test_for_enrichment <- "TaxonomyRank2"
# colname_test_for_enrichment <- "TaxonomyRank3"
# colname_test_for_enrichment <- "TaxonomyRank4"
# colname_test_for_enrichment <- "Neurotransmitter_class"
# 
# # Region
# df.dataset <- df.ldsc_cts.region %>% filter(Class=="Neurons")
# colname_test_for_enrichment <- "Region"
# ### ===== END 'dataset' =====
# 
# ### Set values
# df.dataset <- df.dataset %>% mutate(p.value.mlog10 = -log10(p.value))
# colname_statistic <- "p.value.mlog10"
# # SWITCH: 'Class' or 'category' - or any other column
# # colname_test_for_enrichment <- "Class" 


# ### Run Enrichment
# df.wtest_res <- data.frame(enrichment_val=sapply(df.dataset %>% 
#                                                    filter(!is.na(!!rlang::sym(colname_test_for_enrichment))) %>% # make sure we don't get any NA values
#                                                    distinct(!!rlang::sym(colname_test_for_enrichment)) %>% 
#                                                    pull(), 
#                                                  do_enrichment_test_category, 
#                                                  df=df.dataset, 
#                                                  colname_test_for_enrichment,
#                                                  colname_statistic))
# ### Add number of observations in each group
# df.wtest_res <- df.wtest_res %>% 
#   rownames_to_column(var="category") %>% 
#   as.tibble() %>%
#   arrange(enrichment_val) %>% 
#   left_join(df.dataset %>% count(!!sym(colname_test_for_enrichment)), by=c(category=colname_test_for_enrichment))
# 
# df.wtest_res
# 
