################ SYNOPSIS ###################
# Functions to rename dataset annotation names


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
# ============================== READ METADATA =============================== #
# ======================================================================= #

get_metadata <- function(dataset_prefix) {
  ### DESCRIPTION: function to load expression meta-data
  
  if (dataset_prefix == "campbell_lvl1") {
    file.metadata <- here("/data/expression/hypothalamus_campbell/campbell_lvl1.cell_type_metadata.csv")
    df.metadata <- suppressMessages(read_csv(file.metadata)) %>% rename(annotation = cell_type_all_lvl1)
    df.metadata <- df.metadata %>% mutate(color_by_variable = taxonomy_lvl1)
  } else if (dataset_prefix == "campbell_lvl2") {
    file.metadata <- here("/data/expression/hypothalamus_campbell/campbell_lvl2.cell_type_metadata.csv")
    df.metadata <- suppressMessages(read_csv(file.metadata)) %>% rename(annotation = cell_type_all_lvl2)
    df.metadata <- df.metadata %>% mutate(color_by_variable = taxonomy_lvl1)
  } else if (dataset_prefix == "mousebrain_all") {
    cols_metadata_keep <- c("ClusterName",
                            "Class",
                            "Description",
                            "NCells",
                            "Neurotransmitter",
                            "Probable_location",
                            "Region",
                            "TaxonomyRank1",
                            "TaxonomyRank2",
                            "TaxonomyRank3",
                            "TaxonomyRank4")
    file.metadata <- here("/data/expression/mousebrain/mousebrain-agg_L5.metadata.csv") # PT formatted/cleaned meta-data
    # df.metadata <- read_csv(file.metadata) %>% select(cols_metadata_keep) %>% rename(annotation = ClusterName)
    df.metadata <-  suppressMessages(read_csv(file.metadata))
    df.metadata <- df.metadata %>% mutate(color_by_variable = Class)
  } else if (dataset_prefix == "tabula_muris") {
    file.metadata <- here("/data/expression/tabula_muris/tabula_muris_facs.tissue_celltype.celltype_metadata.csv")
    df.metadata <- suppressMessages(read_csv(file.metadata))
    df.metadata <- df.metadata %>% mutate(annotation = tissue_celltype)
    df.metadata <- df.metadata %>% mutate(color_by_variable = tissue)
    df.metadata <- df.metadata %>% mutate(annotation = stringr::str_replace_all(annotation, pattern="\\s+", replacement="_"))
  } else {
    stop("wrong dataset_prefix")
  }
  return(df.metadata)
}
# ======================================================================= #
# =============================== TM RENAME FUNCTION ================================= #
# ======================================================================= #


utils.rename_annotations.tabula_muris <- function(annotations, style, check_all_matches=F) {
  ### DESCRIPTION
  # Rename Tabula Muris annotation names (for publication)
  ### INPUT
  # annotations:      a vector of annotations. 
  #                   the annotation names should be cleaned without spaces
  #                   e.g. for tabula_muris "Brain_Non-Myeloid.oligodendrocyte_precursor_cell"
  # style:            style.
  # check_all_matches:  if true, the function will raise an exception if not all input annotations were found in the database.
  ### OUTPUT
  # a vector of formatted annotation names. Only matching names will be updated.
  # the length of the output vector is equal to the input vector.
  
  allowed_styles <- c("celltype", "tissue - celltype", "celltype (tissue)")
  if (!style %in% allowed_styles) {
    stop("Got wrong style. Allowed styles are: [%s]", paste(allowed_styles, collapse=", "))
  }
  
  ### Get metadata
  df.metadata <- get_metadata("tabula_muris")
  ### Clean metadata: tissues
  df.metadata <- df.metadata %>% mutate(tissue = case_when(
    tissue == "Brain_Myeloid" ~ "Brain",
    tissue == "Brain_Non-Myeloid" ~ "Brain",
    TRUE ~ as.character(tissue)
    )
  )
  ### Format metadata
  df.metadata <- df.metadata %>% mutate(tissue = stringr::str_replace_all(tissue, pattern="_", replacement=" "))
  df.metadata <- df.metadata %>% mutate(cell_type = stringr::str_replace_all(cell_type, pattern="_", replacement=" "))

  ### Check that all are present in database
  bool.matches <- annotations %in% df.metadata$annotation
  if (check_all_matches && !all(bool.matches)) {
    annotations_not_found <- annotations[!bool.matches]
    stop(sprintf("Some annotations not found: [%s]", paste(annotations_not_found, collapse=",")))
  }
  
  ### Format output conditional on style
  df.metadata.select <- df.metadata %>% filter(annotation %in% annotations)
  if (style == "celltype") {
    df.out_fmt <- df.metadata.select %>% 
      mutate(fmt=cell_type)
  } else if (style=="tissue - celltype") {
    df.out_fmt <- df.metadata.select %>% 
      mutate(fmt=paste(tissue, " - ", cell_type, sep=""))
  } else if (style=="celltype (tissue)") {
    df.out_fmt <- df.metadata.select %>% 
      mutate(fmt=paste(cell_type, "(", tissue, ")", sep=""))
  } else {
    stop("Internal function error")
  }
  
  ### Join renamed with input
  df.all <- tibble(annotation=annotations) # this ensures that the return value has the same length as the input
  df.all <- df.all %>% left_join(df.out_fmt, by="annotation")
  ### Add fmt column for non-matches
  df.all <- df.all %>% mutate(fmt=if_else(is.na(fmt), annotation, fmt)) 
  # ^ non-matches have NA values in the 'fmt' column. we replace na with the input annotation name
  annotations_fmt <- df.all %>% pull(fmt) # returns vector
  return(annotations_fmt)
}

### EXAMPLE CALL
# style <- "tissue - celltype"
# annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5", "Brain_Non-Myeloid.neuron", "Brain_Non-Myeloid.oligodendrocyte_precursor_cell")
# annotations_fmt <- utils.rename_annotations.tabula_muris(annotations, style, check_all_matches=F)
# annotations_fmt


# ======================================================================= #
# ========================= CAMPBELL RENAME FUNCTION ===================== #
# ======================================================================= #


utils.rename_annotations.campbell2015 <- function(annotation_ids, style, check_all_matches=F) {
  stop("NOT IMPLEMENTED YET")
}

# ======================================================================= #
# ========================= MOUSEBRAIN FUNCTION ===================== #
# ======================================================================= #


get_annotations.mousebrain.hypothalamus <- function() {
  print("Will fetch cell-types annotated with 'Hypothalamus' in Region")
  df.metadata <- get_metadata("mousebrain_all")
  
  df.mb_hypo <- df.metadata %>% filter(grepl("Hypothalamus", Region, ignore.case=TRUE))
  ### ALTERNATIVE: Separate rows [also works and gives same result]
  # df.mb_hypo <- df.metadata %>% 
  #   separate_rows(Region, sep = ",") %>% 
  #   mutate(Region = str_trim(Region, side="both")) %>% # trim any leading/trailing whitespace to avoid e.g. " Striatum ventral" and "Striatum ventral" being separate Regions.
  #   filter(Region=="Hypothalamus")
  #
  stopifnot(length(unique(df.mb_hypo$annotation))==nrow(df.mb_hypo)) # ensure all annotations are unique
  annotations.hypo <- df.mb_hypo %>% pull(annotation)
  print(sprintf("Returning vector of n=%s annotations", length(annotations.hypo)))
  return(annotations.hypo)
}

