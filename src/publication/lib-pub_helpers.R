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


# ======================================================================= #
# ================================= MISC ================================ #
# ======================================================================= #

get_scrna_seq_dataset_prefixes <- function(selection=c("all", "brain", "hypo")) {
  selection <- match.arg(selection)
  ### These names should match the prefixes in /out/es/*es_obj.RData for the "final publication" datasets
  data_prefixes <- c("campbell2017_lvl2",
                     "chen2017",
                     "kimVMH2019_smartseq",
                     "mikkelsen2019",
                     "moffitt2018",
                     "mousebrain",
                     "romanov2017",
                     "tabula_muris")
  if (selection == "all") {
    # do nothing
  } else if (selection =="brain") {
    data_prefixes <- data_prefixes[data_prefixes!="tabula_muris"]
  } else if (selection =="hypo") {
    data_prefixes <- data_prefixes[!data_prefixes %in% c("tabula_muris","mousebrain")]
  }
  message(sprintf("selection = %s; returning n=%s data prefixes", selection, length(data_prefixes)))
  message(paste(data_prefixes, collapse=" | "))
  return(data_prefixes)
}

# ======================================================================= #
# ========================= PRIORITIZED CELL-TYPES ====================== #
# ======================================================================= #

get_prioritized_annotations_bmi <- function(dataset) {
  ### BMI_UKBB_Loh2018 dataset FDR significant
  allowed.dataset <- c("mousebrain", "tabula_muris", "campbell2017_lvl2")
  if (!dataset %in% allowed.dataset) {
    stop("Got wrong dataset argument. Allowed values are: [%s]", paste(allowed.dataset, collapse=", "))
  }
  if (dataset == "mousebrain") {
    # annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
    # annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","DEGLU4","TEINH12") # No TEGLU4
    annotations <- c("DEINH3","MEGLU11","MEINH2","MEGLU10","DEGLU5","MEGLU1","MEGLU3","TEGLU23","DEGLU4","TEGLU19","HBSER5","HBGLU2","MEGLU2","TEGLU17","MEGLU7","MEINH3","HBGLU5","TEINH12","TEINH1","TEGLU4","TEGLU21","HBINH8")
  } else if (dataset == "tabula_muris") {
    annotations <- c("Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte_precursor_cell")
  } else if (dataset == "campbell2017_lvl2") {
    annotations <- c("n29.Nr5a1-Adcyap1") # Cleaned name is n29.Nr5a1/Bdnf
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
  
  allowed.dataset <- c("mousebrain", "tabula_muris", "campbell2017_lvl2")
  if (!dataset %in% allowed.dataset) {
    stop("Got wrong dataset argument. Allowed values are: [%s]", paste(allowed.dataset, collapse=", "))
  }
  if (dataset == "mousebrain") {
    # c(name=value)
    ### thesis version
    # annotations <- c("TEGLU23"="#E31A1C", # HC
    #                  # "TEGLU4"="#FF7F00", # CORTEX
    #                  "TEGLU17"="#FF7F00", # CORTEX
    #                  "TEINH12"="#B15928", # HC/CORTEX
    #                  "DEGLU5"="#1F78B4", # MIDBRAIN
    #                  "DEINH3"="#33A02C", # HYPOTHALAMUS
    #                  "MEGLU1"="#1F78B4", # MIDBRAIN
    #                  "MEINH2"="#1F78B4", # MIDBRAIN
    #                  "MEGLU11"="#1F78B4", # MIDBRAIN
    #                  "MEGLU10"="#1F78B4", # MIDBRAIN
    #                  "DEGLU4"="#6A3D9A" # THALAMUS
    #                  )
    ### v1 (n=22)
    # annotations <- c("TEGLU23"="#E31A1C", # HC
    #                  "TEGLU4"="#FF7F00", # CORTEX
    #                  "TEGLU17"="#FF7F00", # CORTEX
    #                  "TEGLU19"="#FF7F00", # CORTEX
    #                  "TEGLU21"="#E31A1C", # HC
    #                  "TEINH12"="#B15928", # HC/CORTEX
    #                  "HBINH8"="#5ab4ac", # MEDULLA
    #                  "DEGLU4"="#6A3D9A", # THALAMUS
    #                  "DEINH3"="#33A02C", # HYPOTHALAMUS
    #                  "DEGLU5"="#1F78B4", # MIDBRAIN
    #                  "MEGLU1"="#1F78B4", # MIDBRAIN
    #                  "MEINH2"="#1F78B4", # MIDBRAIN
    #                  "MEGLU11"="#1F78B4", # MIDBRAIN
    #                  "MEGLU10"="#1F78B4", # MIDBRAIN
    #                  "MEGLU2"="#1F78B4", # MIDBRAIN
    #                  "MEGLU3"="#1F78B4", # MIDBRAIN
    #                  "MEGLU7"="#1F78B4", # MIDBRAIN
    #                  "MEINH3"="#1F78B4", # MIDBRAIN
    #                  "TEINH1"="black", # Pallidum
    #                  "HBGLU2"="#5ab4ac", # MEDULLA
    #                  "HBGLU5"="#c51b7d", # PONS
    #                  "HBSER5"="#c51b7d" # PONS
    # ) 
    ### v2 (n=22) | Hindbrain
    annotations <- c("TEGLU23"="#E31A1C", # HC
                     "TEGLU4"="#FF7F00", # CORTEX
                     "TEGLU17"="#FF7F00", # CORTEX
                     "TEGLU19"="#FF7F00", # CORTEX
                     "TEGLU21"="#E31A1C", # HC
                     "TEINH12"="#c51b7d", # HC/CORTEX
                     "HBINH8"="#6A3D9A", # MEDULLA
                     "DEGLU4"="#01665e", # THALAMUS
                     "DEINH3"="#33A02C", # HYPOTHALAMUS
                     "DEGLU5"="#1F78B4", # MIDBRAIN
                     "MEGLU1"="#1F78B4", # MIDBRAIN
                     "MEINH2"="#1F78B4", # MIDBRAIN
                     "MEGLU11"="#1F78B4", # MIDBRAIN
                     "MEGLU10"="#1F78B4", # MIDBRAIN
                     "MEGLU2"="#1F78B4", # MIDBRAIN
                     "MEGLU3"="#1F78B4", # MIDBRAIN
                     "MEGLU7"="#1F78B4", # MIDBRAIN
                     "MEINH3"="#1F78B4", # MIDBRAIN
                     "TEINH1"="#b35806", # Pallidum
                     "HBGLU2"="#6A3D9A", # MEDULLA
                     "HBGLU5"="#6A3D9A", # PONS
                     "HBSER5"="#6A3D9A" # PONS
    ) 
  } else if (dataset == "tabula_muris") {
    # #FF0000 == red
    annotations <- c("Brain_Non-Myeloid.neuron"="#FF0000", 
                     "Brain_Non-Myeloid.oligodendrocyte_precursor_cell"="#FF0000")
  } else if (dataset == "hypothalamus") {
    annotations <- c("n29.Nr5a1-Adcyap1"="#FF0000")
    # i9_Gaba
    # n29.Nr5a1-Adcyap1
    # Glut_5
    # e8_Cck_Ebf3
  } else if (dataset == "campbell2017_lvl2") {
    stop("Campbell option no longer supported")
    annotations <- c("n29.Nr5a1-Adcyap1"="#FF0000") # Cleaned name is n29.Nr5a1/Bdnf
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
  # LEGEND WILL BE ORDERED ACCORDING TO THIS ORDER
  colormap.region <- c("Cortex"="#FF7F00",
                     "Hippocampus/Cortex"="#c51b7d", # prev. brown==> B15928
                     "Hippocampus"= "#E31A1C", 
                     "Hindbrain"="#6A3D9A",
                     "Midbrain"="#1F78B4",
                     "Pallidum"="#b35806", # black
                     "Thalamus"="#01665e",
                     "Subthalamus"="#33A02C" # "Hypothalamus"="#33A02C",
                     # "Medulla"="#5ab4ac",
                     # "Pons"="#c51b7d",
                     )
  return(colormap.region)
}

