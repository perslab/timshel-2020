############### SYNOPSIS ###################
# Helper functions for publication figure related to colors and annotation selections


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


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
  allowed.dataset <- c("mousebrain", "tabula_muris", "hypothalamus")
  if (!dataset %in% allowed.dataset) {
    stop("Got wrong dataset argument. Allowed values are: [%s]", paste(allowed.dataset, collapse=", "))
  }
  if (dataset == "mousebrain") {
    # annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
    # annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","DEGLU4","TEINH12") # No TEGLU4
    annotations <- c("DEINH3","MEGLU11","MEINH2","MEGLU10","DEGLU5","MEGLU1","MEGLU3","TEGLU23","DEGLU4","TEGLU19","HBSER5","HBGLU2","MEGLU2","TEGLU17","MEGLU7","MEINH3","HBGLU5","TEINH12","TEINH1","TEGLU4","TEGLU21","HBINH8")
  } else if (dataset == "tabula_muris") {
    annotations <- c("Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte_precursor_cell")
  } else if (dataset == "hypothalamus") {
    annotations <- c("POA-NEURO66","ARCME-NEURO29","LHA-NEURO20","POA-NEURO21")
    # } else if (dataset == "campbell2017_lvl2") {
    #   stop("Campbell not supported any more")
    #   annotations <- c("n29.Nr5a1-Adcyap1") # Cleaned name is n29.Nr5a1/Bdnf
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
  
  allowed.dataset <- c("mousebrain", "tabula_muris", "hypothalamus")
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
    annotations.raw <- c("POA-NEURO66","ARCME-NEURO29","LHA-NEURO20","POA-NEURO21")
    annotations <- rep("#FF0000", times=length(annotations.raw))
    names(annotations) <- annotations.raw
    # "i9_Gaba"
    # "n29.Nr5a1-Adcyap1"
    # "Glut_5"
    # "e8_Cck_Ebf3"
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
# ============================ MOUSEBRAIN region color ======================= #
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
