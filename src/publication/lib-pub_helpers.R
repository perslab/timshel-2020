############### SYNOPSIS ###################
# Helper functions for publication figures


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ========================= LIBRARY DENPENDENCIES ======================= #
# ======================================================================= #

### These variables are assumed to be available in the global environment
# df.ldsc_cts

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)


# ======================================================================= #
# ================================== UTILS ============================== #
# ======================================================================= #

reorder_factor_levels_to_first_last_pairs <- function(factor_levels) {
  ### reorder vector. works for input vectors of both equal and unequal lengths
  ### this function is useful for adding bigger contrast between neighboring colored categories.
  ### input
  # factor_levels: ordered factor levels (unique). The factor levels should be ordered in the same order as your x/y annotations are ordered
  ### output
  # a charactor vector with the same length and elements as the input factor_levels. 
  # the elements are re-ordered to form 'first-last' pairs
  ### usage
  # df <- df %>% mutate(x = factor(x, levels=reorder_factor_levels_to_first_last_pairs(unique(x)), ordered=T))
  ### theory/explanation for why 'correction' for an input with equal number of elements ('equal length') is needed.
  # 'selected' indexes are marked in '[ ]'
  # equal length without correction: ---> 2 and 4 skipped
  # [1] 2 [3] 4
  #  4 [3] 2 [1]
  # equal length WITH correction: ---> ok, works
  # [1] 2 [3] 4
  #  X [4] 3 [2] 1 --> correction: revered indexes are shifted one place
  # unequal length without correction: --> ok, no correction needed
  # [1] 2 [3] 4  [5]
  #  5 [4] 3 [2] 1
  if (any(duplicated(factor_levels))) {
    stop("Duplicated values in input factor_levels")
  }
  if (length(factor_levels) %% 2 == 0) { # equal length input
    # print("Running with correction mode")
    flag_correction <- TRUE
  } else {
    flag_correction <- FALSE
  }
  idx_levels <- seq(1, length(factor_levels))
  factor_levels.reorder <- vector(mode="character", length=length(factor_levels))
  for (i in idx_levels) {
    if (i %% 2 == 0) {
      if (flag_correction) {
        idx <- c(999, rev(idx_levels))[i] # 999 is a dummy number.
      } else {
        idx <- rev(idx_levels)[i]
      }
      # print(sprintf("rev idx = %s", idx))
    } else {
      idx <- idx_levels[i]
      # print(sprintf("fwd idx = %s", idx))
    }
    factor_levels.reorder[i] <- factor_levels[idx]
    # print(sprintf("factor_levels[idx]=%s", factor_levels[idx]))
  }
  ### Safety checks (unit testing)
  stopifnot(length(factor_levels.reorder) == length(factor_levels))
  if (any(duplicated(factor_levels.reorder))) {
    stop("Duplicated values in output")
  }
  return(factor_levels.reorder)
}

# ======================================================================= #
# ============================= LOADER FUNCTIONS ======================== #
# ======================================================================= #

format_cellect_ldsc_results <- function(df.ldsc_cts) {
  ### Rename 
  df.ldsc_cts <- df.ldsc_cts %>% rename(annotation=annotation,
                            specificity_id=specificity_id,
                            estimate=beta,
                            std.error=beta_se,
                            p.value=pvalue)
  # df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n_obs_es, true=T, false=F),
  #                                       p.value.adj = p.value*n_obs_es)
  return(df.ldsc_cts)
}


