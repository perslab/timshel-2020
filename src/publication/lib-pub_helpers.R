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
  ### this function is usefull for adding bigger contrast between neighboring colored categories.
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
# ========================= PRIORITIZED CELL-TYPES ====================== #
# ======================================================================= #

get_prioritized_annotations_bmi <- function(dataset) {
  ### BMI_UKBB_Loh2018 dataset FDR significant
  allowed.dataset <- c("mousebrain", "tabula_muris", "campbell")
  if (!dataset %in% allowed.dataset) {
    stop("Got wrong dataset argument. Allowed values are: [%s]", paste(allowed.dataset, collapse=", "))
  }
  if (dataset == "mousebrain") {
    annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
  } else if (dataset == "tabula_muris") {
    annotations <- c("Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte_precursor_cell")
  } else if (dataset == "campbell") {
    stop("Not yet implemented")
  } else {
    stop("Internal error")
  }
  print(sprintf("Returning n=%s annotation ids", length(annotations)))
  return(annotations)
}

# ======================================================================= #
# ============================ CELL-TYPES: COLORS ======================= #
# ======================================================================= #


get_color_mapping.prioritized_annotations_bmi <- function(dataset) {
  ### EXAMPLE USAGE
  # colormap.annotations <- get_prioritized_annotations_color_mapping(dataset="mousebrain")
  # ggplot(...) + scale_color_manual(values=colormap.annotation) # a set of aesthetic values to map data values to. the values will be matched based on the names.
  
  allowed.dataset <- c("mousebrain", "tabula_muris", "campbell")
  if (!dataset %in% allowed.dataset) {
    stop("Got wrong dataset argument. Allowed values are: [%s]", paste(allowed.dataset, collapse=", "))
  }
  if (dataset == "mousebrain") {
    # "Cortex"="#FF7F00",
    # "Hippocampus/Cortex"="#B15928",
    # "Hippocampus"="#E31A1C",
    # "Thalamus"="#6A3D9A",
    # "Midbrain"="#1F78B4",
    # "Hypothalamus"="#33A02C"
    ### ORDERED BY HIERARCHICAL CLUSTERING OF ESmu
    annotations <- c("TEGLU23"="#E31A1C", # HC
                     "TEGLU4"="#FF7F00", # CORTEX
                     "TEGLU17"="#FF7F00", # CORTEX
                     "TEINH12"="#B15928", # HC/CORTEX
                     "DEGLU5"="#1F78B4", # MIDBRAIN
                     "DEINH3"="#33A02C", # HYPOTHALAMUS
                     "MEGLU1"="#1F78B4", # MIDBRAIN
                     "MEINH2"="#1F78B4", # MIDBRAIN
                     "MEGLU11"="#1F78B4", # MIDBRAIN
                     "MEGLU10"="#1F78B4", # MIDBRAIN
                     "DEGLU4"="#6A3D9A" # THALAMUS
                     ) # c(name=value)
  } else if (dataset == "tabula_muris") {
    # #FF0000 == red
    annotations <- c("Brain_Non-Myeloid.neuron"="#FF0000", 
                     "Brain_Non-Myeloid.oligodendrocyte_precursor_cell"="#FF0000")
  } else if (dataset == "campbell") {
    stop("Not yet implemented")
  } else {
    stop("Internal error")
  }
  print(sprintf("Returning n=%s annotation ids. x = colors; names(x) = annotations", length(annotations)))
  return(annotations)
}


# ======================================================================= #
# ============================ MOUSEBRAIN FUNCTIONS ======================= #
# ======================================================================= #


get_color_mapping.mb.region <- function() {
  colormap.region <- c("Cortex"="#FF7F00",
                     "Hippocampus/Cortex"="#B15928",
                     "Hippocampus"="#E31A1C",
                     "Thalamus"="#6A3D9A",
                     "Midbrain"="#1F78B4",
                     "Hypothalamus"="#33A02C")
  return(colormap.region)
}

