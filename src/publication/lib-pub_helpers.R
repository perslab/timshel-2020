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

get_prioritized_annotations_color_mapping <- function(dataset) {
  ### EXAMPLE USAGE
  # colormap.annotations <- get_annotation_color_mapping(dataset="mousebrain")
  # ggplot(...) + scale_color_manual(values=colormap.annotation)
  
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
    annotations <- c("TEGLU4"="#FF7F00", # CORTEX
                     "TEGLU17"="#FF7F00", # CORTEX
                     "TEGLU23"="#E31A1C", # HC
                     "TEINH12"="#B15928", # HC/CORTEX
                     "DEGLU4"="#6A3D9A", # THALAMUS
                     "DEINH3"="#33A02C", # HYPOTHALAMUS
                     # MIDBRAIN
                     "MEGLU1"="#1F78B4",
                     "MEINH2"="#1F78B4",
                     "DEGLU5"="#1F78B4",
                     "MEGLU10"="#1F78B4",
                     "MEGLU11"="#1F78B4"
                     )
  } else if (dataset == "tabula_muris") {
    stop("Not yet implemented")
  } else if (dataset == "campbell") {
    stop("Not yet implemented")
  } else {
    stop("Internal error")
  }
  print(sprintf("Returning n=%s annotation ids. x = annotations; names(x) = colors.", length(annotations)))
  return(annotations)
}


