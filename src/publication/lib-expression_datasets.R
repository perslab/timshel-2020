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
# ============================== READ METADATA =============================== #
# ======================================================================= #

get_metadata <- function(dataset_prefix) {
  ### DESCRIPTION: function to load expression meta-data
  
  if (dataset_prefix == "campbell2017_lvl1") {
    # stop(sprintf("dataset_prefix '%s' is no longer supported (December 2019).", dataset_prefix))
    file.metadata <- here("/data/expression/campbell2017/campbell2017_lvl1.cell_type_metadata.csv")
    df.metadata <- suppressMessages(read_csv(file.metadata)) %>% rename(annotation = cell_type_all_lvl1)
    df.metadata <- df.metadata %>% mutate(color_by_variable = taxonomy_lvl1)
  } else if (dataset_prefix == "campbell2017_lvl2") {
    # stop(sprintf("dataset_prefix '%s' is no longer supported (December 2019).", dataset_prefix))
    file.metadata <- here("/data/expression/campbell2017/campbell2017_lvl2.cell_type_metadata.csv")
    df.metadata <- suppressMessages(read_csv(file.metadata)) %>% rename(annotation = cell_type_all_lvl2)
    df.metadata <- df.metadata %>% mutate(color_by_variable = taxonomy_lvl1)
  } else if (dataset_prefix == "hypothalamus") {
    file.metadata <- here("data/expression/hypothalamus/hypothalamus_metadata.csv")
    df.metadata <- suppressMessages(read_csv(file.metadata))
    file.tax <- here("data/expression/hypothalamus/hypothalamus_taxonomy.csv")
    df.tax <- suppressMessages(read_csv(file.tax))
    df.metadata <- df.metadata %>% left_join(df.tax, by="taxonomy_lvl2")
    # df.metadata <- df.metadata %>% select(annotation, specificity_id, annotation_prefix, annotation_fmt, annotation_fmt_prefix)
    df.metadata <- df.metadata %>% mutate(color_by_variable = annotation_prefix)
  } else if (dataset_prefix == "mousebrain") {
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
    file.metadata <- here("/data/expression/mousebrain/mousebrain.metadata.csv") # PT formatted/cleaned meta-data
    # df.metadata <- read_csv(file.metadata) %>% select(cols_metadata_keep) %>% rename(annotation = ClusterName)
    df.metadata <-  suppressMessages(read_csv(file.metadata))
    df.metadata <- df.metadata %>% mutate(color_by_variable = Class) # not needed anymore?
    ### Add Region_fmt [specifically made for prioritized cell-types]
    df.metadata <- df.metadata %>% mutate(Region_fmt = case_when(
      annotation == "DEINH3" ~ "Subthalamus", # Subthalamic nucleus
      Region %in% c("Pons", "Medulla") ~ "Hindbrain",
      Region == "Midbrain dorsal" ~ "Midbrain",
      Region == "Midbrain dorsal,Midbrain ventral" ~ "Midbrain",
      Region == "Hippocampus,Cortex" ~ "Hippocampus/Cortex",
      TRUE ~ as.character(Region))
    )
    
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
# ========================= HYPOTHALAMUS RENAME FUNCTION ===================== #
# ======================================================================= #

utils.rename_annotations.hypothalamus <- function(annotations, specificity_ids, check_all_matches=F) {
  #' Function to rename hypothalamus dataset
  #' 
  #' This function is complicated by duplicated annotation names (`GABA1`, `GABA_10`, `GABA_8`, `GABA_9`)
  #' That is why both a vector of annotations and specificity_ids must be given.
  #' SEE utils.rename_annotations.tabula_muris() for other details

  ### TMP debug
  # specificity_ids <- df$specificity_id
  # annotations <- df$annotation
  # check_all_matches <- TRUE
  
  ### Load hypothalamus metadata
  df.metadata <- get_metadata("hypothalamus")
  df.metadata <- df.metadata %>% mutate(annotation_uniq = paste0(specificity_id, "__", annotation)) # needed because some hyp annotations are duplicated across datasets
  
  ### Make unique names
  annotation_uniq = paste0(specificity_ids, "__", annotations)

  ### Check that all are present in database
  bool.matches <- annotation_uniq %in% df.metadata$annotation_uniq
  if (check_all_matches && !all(bool.matches)) {
    annotations_not_found <- annotation_uniq[!bool.matches]
    stop(sprintf("Some annotation_uniq not found: [%s]", paste(annotations_not_found, collapse=",")))
  }
  
  ### Format output [legacy from TM function]
  df.out_fmt <- df.metadata %>% filter(annotation_uniq %in% annotation_uniq) %>% mutate(fmt=annotation_fmt) # annotation_fmt is the unique and formatted name in hyp meta-data
  
  ### Join renamed with input
  df.all <- tibble(annotation_uniq=annotation_uniq, annotation=annotations) # this ensures that the return value has the same length as the input
  df.all <- df.all %>% left_join(df.out_fmt %>% select(-annotation), by="annotation_uniq") # select(-annotation) to avoid annotation.x and annotation.y
  ### Add fmt column for non-matches
  df.all <- df.all %>% mutate(fmt=if_else(is.na(fmt), annotation, fmt)) # *OBS*: notice we use annotation and not annotation_uniq to ensure non-mapped annotations does not get changed in the output
  # ^ non-matches have NA values in the 'fmt' column. we replace na with the input annotation name
  annotations_fmt <- df.all %>% pull(fmt) # returns vector
  return(annotations_fmt)
}

### OLD MANUAL METHOD - WORKS [OK DELETE]
# ### Load hypothalamus metadata with 'clean names'
# df.metadata.hyp <- get_metadata("hypothalamus")
# df.metadata.hyp <- df.metadata.hyp %>% mutate(annotation_uniq = paste0(specificity_id, "__", annotation)) # needed because some hyp annotations are duplicated across datasets
# 
# ### Map hypothalamus annotations to 'clean names'
# df <- df %>% mutate(annotation_uniq = paste0(specificity_id, "__", annotation))
# df <- df %>% mutate(annotation_fmt = df.metadata.hyp$annotation_fmt[match(annotation_uniq, df.metadata.hyp$annotation_uniq)]) # mousebrain annotation will get NA values
# df <- df %>% mutate(annotation_fmt = if_else(is.na(annotation_fmt), annotation, annotation_fmt)) # replace NA with original value (cell_type)


# ======================================================================= #
# ========================= MOUSEBRAIN HYP FUNCTION ===================== #
# ======================================================================= #


get_annotations.mousebrain.hypothalamus <- function() {
  print("Will fetch cell-types annotated with 'Hypothalamus' in Region")
  df.metadata <- get_metadata("mousebrain")
  
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



# ======================================================================= #
# ================ HYPOTHALAMUS: load combined ESmu matrix ============== #
# ======================================================================= #
### OUTPUT: df.es
# A tibble: 6 x 348
# gene  `ARCME-NEURO1` `ARCME-NEURO2` `ARCME-NEURO3` `ARCME-NEURO4` `ARCME-NEURO5` `ARCME-NEURO6` `ARCME-NEURO7` `ARCME-NEURO8`
# ENSG…          0              0             0              0.0306          0.157         0             0.119            0    
# ENSG…          0              0             0              0               0             0             0.185            0    
# ENSG…          0.155          0.174         0.0614         0.187           0             0.0702        0.200   

get_combined_hypo_es <- function(merge=c("inner", "full")) {
  ### Meta-data
  df.metadata <- get_metadata("hypothalamus")
  df.metadata <- df.metadata %>% mutate(annotation_uniq = paste0(specificity_id, "__", annotation)) # needed because some hyp annotations are duplicated across datasets
  
  ### ESmu matrix
  file.es <- here("out/es", paste0(get_scrna_seq_dataset_prefixes("hypo"), ".mu.csv.gz"))
  list.df.es <- map(file.es, read_csv) # genes x cell-types
  ### Rename ESmu annotations to avoid problem with duplicates during merge 
  names(list.df.es) <- get_scrna_seq_dataset_prefixes("hypo") # names is garantueed to be in same order as above
  list.df.es.prefixed <- list.df.es %>% imap(.f=function(df,specificity_id){colnames(df)[-1] <- paste0(specificity_id, "__", colnames(df)[-1]); df}) # genes x cell-types.
  # ^ x[-1]: first column is "gene"
  
  ### Merge
  df.es <- plyr::join_all(list.df.es.prefixed, by="gene", type=merge, match="first") %>% as.tibble() # REF: https://stackoverflow.com/a/32066419/6639640. match argument should not matter (only speed)
  # plyr::join_all type: left, right, inner or full
  # *inner*: only overlapping genes will be used for clustering
  # *full*: used for exporting combined ESmu matrix
  
  ### Map hypothalamus annotations to annotation_fmt
  annotation_uniq <- colnames(df.es)[-1] # exclude 'gene' column
  annotation_fmt <- df.metadata$annotation_fmt[match(annotation_uniq, df.metadata$annotation_uniq)]
  # df.tmp <- tibble(x=annotation_fmt, y=annotation_uniq) # test to show that the matching works
  colnames(df.es)[-1] <- annotation_fmt # ALTERNATIVE: use deframe() and rename(!!!x)
  return(df.es)
}
