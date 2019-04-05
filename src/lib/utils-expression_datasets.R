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
# =============================== TM RENAME FUNCTION ================================= #
# ======================================================================= #

# ### tissue: rename + remove some text because they contain too few cell-types
# df.plot.tax_order %>% count(tissue, sort=T) # ---> potentially filter : n_annotations_in_tax >= 5
# df.plot.tax_order <- df.plot.tax_order %>% mutate(tissue = case_when(
#   tissue == "Brain_Myeloid" ~ "Brain",
#   tissue == "Brain_Non-Myeloid" ~ "Brain",
#   TRUE ~ stringr::str_replace(tissue, pattern="_", replacement=" ") # TRUE ~ as.character(tissue)
# )
# )


# # Self-defined label formatting function
# tm_annotation_formatter <- function(annotation_tissue_celltype) {
#   # INPUT: annotation_tissue_celltype, string e.g. 'Marrow.Slamf1-negative_multipotent_progenitor_cell'
#   # USAGE 1: df %>% mutate
#   # USAGE 2: ggplot() + scale_y_discrete(label=tm_annotation_formatter)
#   label <- stringr::str_split_fixed(annotation_tissue_celltype, pattern="\\.", n=Inf)[,2]
#   label <- stringr::str_replace_all(label, pattern="-", " - ") # add space to any hyphens
#   label <- stringr::str_replace_all(label, pattern="_", " ") # convert _ to space
#   # label <- stringr::str_to_sentence(label) # title case
#   substr(label, 1, 1) <- toupper(substr(label, 1, 1)) # capitalize first character
#   return(label)
# }


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
    TRUE ~ stringr::str_replace(tissue, pattern="_", replacement=" ") # TRUE ~ as.character(tissue)
    )
  )
  ### Note: "cell_type" column is already clean: Bergmann glial cell", "brain pericyte"

  
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
# ============================== READ METADATA =============================== #
# ======================================================================= #

get_metadata <- function(dataset_prefix) {
  ### DESCRIPTION: function to load expression meta-data
  
  if (dataset_prefix == "campbell_lvl1") {
    file.metadata <- here("/data/expression/hypothalamus_campbell/campbell_lvl1.cell_type_metadata.csv")
    df.metadata <- read_csv(file.metadata) %>% rename(annotation = cell_type_all_lvl1)
    df.metadata <- df.metadata %>% mutate(color_by_variable = taxonomy_lvl1)
  } else if (dataset_prefix == "campbell_lvl2") {
    file.metadata <- here("/data/expression/hypothalamus_campbell/campbell_lvl2.cell_type_metadata.csv")
    df.metadata <- read_csv(file.metadata) %>% rename(annotation = cell_type_all_lvl2)
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
    df.metadata <- read_csv(file.metadata)
    df.metadata <- df.metadata %>% mutate(color_by_variable = Class)
  } else if (dataset_prefix == "tabula_muris") {
    file.metadata <- here("/data/expression/tabula_muris/tabula_muris_facs.tissue_celltype.celltype_metadata.csv")
    df.metadata <- read_csv(file.metadata)
    df.metadata <- df.metadata %>% mutate(annotation = tissue_celltype)
    df.metadata <- df.metadata %>% mutate(color_by_variable = tissue)
    df.metadata <- df.metadata %>% mutate(annotation = stringr::str_replace_all(annotation, pattern="\\s+", replacement="_"))
  } else {
    stop("wrong dataset_prefix")
  }
  return(df.metadata)
}

