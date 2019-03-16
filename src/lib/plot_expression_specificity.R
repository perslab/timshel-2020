

library(tidyverse)
library(here)

library(gghighlight)
library(ggrepel)
library(viridis)

# ===================================================================================================== #
# ================================ Hybrid: annotation or genes centric ================================ #
# ===================================================================================================== #

.get_es_data <- function(sem_obj, annotations=NULL, genes=NULL) {
  ### General purpose ES data extractor.
  ### *** FUNCTION NOT FINISHED - it becomes too complicated/inefficient. ***
  
  ### DEVELOPMENT
  # sem_obj
  # annotations <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12")
  # genes <- c("AGRP", "NPY")
  
  if ( is.null(annotations) & is.null(genes) ) {
    stop("You must provide one of 'annotations' or 'genes' arguments")
  } 
  # if ( !is.null(annotations) & !is.null(genes) ) {
  #   stop("You must provide ONLY one of 'annotations' or 'genes' arguments")
  # }
  # 
  if (!is.null(annotations) & is.null(sem_obj[["group_by_annotation.sem"]][[annotation]])) {
    stop(sprintf("Annotation = '%s' not found in data", annotation))
  }
  
  ### Get 'default' values
  if (is.null(annotations)) {
    annotations <- sem_obj$annotations
  }
  if (is.null(genes)) {
    genes <- sem_obj$genes
  }
  
  # bind_rows(sem_obj[["group_by_annotation.sem"]][annotations], .id="annotation")
}


# ==================================================================================== #
# ================================ Annotation centric ================================ #
# ==================================================================================== #
get_es_annotation_centric <- function(sem_obj, annotation) {
  ### OBS: this function only supports extracting one annotation.
  
  if (is.null(sem_obj[["group_by_annotation.sem"]][[annotation]])) {
    stop(sprintf("Annotation = '%s' not found in data", annotation))
  }
  ### Get data
  df.es <- sem_obj[["group_by_annotation.sem"]][[annotation]] %>% 
    mutate(gene=sem_obj$genes) %>% 
    hs_add_gene_symbol_from_ensembl_ids(colname_geneids_from="gene", colname_geneids_to="gene_name") # map gene
  
  ### Add SEM mean
  df.es <- df.es %>% 
    mutate(
      `mean nES` = sem_obj$sem_meta$mean[[annotation]] # returns a (unnamed) numeric vector. This works even though sem_obj$sem_meta$mean is a tibble
    )
  
  ### Wrangle data
  df.es.mod <- df.es %>% 
    mutate(gene_idx_fixed = seq_len(n())) %>%
    gather(key="es_metric", value="es_weight", -gene, -gene_name, -gene_idx_fixed) %>%
    group_by(es_metric) %>%
    mutate(gene_idx_sorted_within_metric = rank(es_weight, ties.method = "first")) %>%
    ungroup()
  
  return(df.es.mod)
  #   gene            gene_name gene_idx_fixed es_metric es_weight gene_idx_sorted_within_metric
  #   <chr>           <chr>              <int> <chr>         <dbl>                         <int>
  #   ENSG00000141668 CBLN2                  1 DE t-test    -0.885                          1040
  #   ENSG00000204624 DISP3                  2 DE t-test    -0.653                          7517
  #   ENSG00000187848 P2RX2                  3 DE t-test    -1.33                            211
  #   ENSG00000171522 PTGER4                 4 DE t-test    -0.722                          4716
}


es_plot_annotation_centric <- function(df.es, genes_highlight) {
  # df.es <- df.es %>% mutate(gene_label = paste0(gene_name, "(", round(es_weight,2), ")") ) # looks too complication
  # TODO: try to add geom_density2d() to the plot
  
  ### rename
  df.es <- df.es %>% mutate(es_metric = case_when(
    es_metric == "ges" ~ "GES",
    es_metric == "specificity" ~ "Specificity",
    es_metric == "si" ~ "SI",
    es_metric == "tstat" ~ "DE t-test",
    TRUE ~ as.character(es_metric))
  )
  
  df.es <- df.es %>% mutate(gene_label = gene_name)
  p <- ggplot(df.es, aes(x=gene_idx_fixed, y=es_weight)) + 
    geom_point(color="red") + 
    # geom_density2d() +
    gghighlight(gene_name %in% genes_highlight, label_key=gene_label) + 
    facet_wrap(~es_metric, scales = "free") + 
    labs(x="Gene", y="ES weight") + 
    # theme_light() # this looks ok as well.
    theme_classic()
  return(p)
}

# ============================================================================== #
# ================================ Gene centric ================================ #
# ============================================================================== #

get_es_gene_centric <- function(sem_obj, genes_select) {
  # this function supports multiple genes (genes_select as vector)
  
  ### DEVELOPMENT
  # sem_obj
  # genes_select <- c("AGRP", "NPY")
  
  ### Get index of genes
  df.genes <- tibble(gene=sem_obj$genes) %>% hs_add_gene_symbol_from_ensembl_ids(colname_geneids_from="gene", colname_geneids_to="gene_name")
  genes_match.idx <- match(genes_select, df.genes$gene_name)
  if (anyNA(genes_match.idx)) {
    print("Warning: one or more genes in genes_select could not be found.")
    genes_not_found <- genes_select[is.na(genes_match.idx)]
    print(sprintf("The following genes could not be found: %s", paste0(genes_not_found, collapse=";")))
    print("These genes will be excluded from return data frame.")
    genes_match.idx <- genes_match.idx[!is.na(genes_match.idx)] # exclude NAs from match
  }
  genes_match.names <- df.genes$gene_name[genes_match.idx] # extract after excluding NAs
  
  list.es <- list()
  for (annotation in sem_obj$annotations) {
    df.es <- sem_obj[["group_by_annotation.sem"]][[annotation]] %>% slice(genes_match.idx)
    es_mean <- sem_obj$sem_meta$mean %>% slice(genes_match.idx) %>% pull(!!sym(annotation)) # vector
    list.es[[annotation]] <- df.es %>% mutate(
      es_mean=es_mean,
      gene_name=genes_match.names)
  }
  df.es <- bind_rows(list.es, .id="annotation")
  # annotation  tstat   ges        si specificity es_mean gene_nanme 
  # <chr>       <dbl> <dbl>     <dbl>       <dbl>   <dbl> <chr>
  # ABC        -0.869 0     0.0000664    0         0      AGRP 
  # ABC        -0.966 0.160 0.660        0.000520  0      NPY  
  
  df.es.gather <- df.es %>% gather(key="es_metric", value="es_weight", -gene_name, -annotation)
  df.es.gather <- df.es.gather %>%
    group_by(gene_name, es_metric) %>%
    mutate(gene_idx_sorted_within_metric = rank(es_weight, ties.method = "first")) %>%
    # mutate(gene_idx_sorted_within_metric = rank(es_weight)) %>%
    ungroup()
  
  
  return(df.es.gather)
  # annotation gene_name es_metric es_weight
  # <chr>      <chr>     <chr>         <dbl>
  # ABC        AGRP      tstat        -0.869
  # ABC        NPY       tstat        -0.966
  # ACBG       AGRP      tstat        -0.869
}


get_mean_nes_gene_centric <- function(sem_obj, genes_select) {
  ### OBS: this function uses only sem_meta$mean data
  # this function supports multiple genes (genes_select as vector)
  
  ### Get data
  df.es_mean <- sem_obj$sem_meta$mean %>% 
    mutate(gene=sem_obj$genes) %>% 
    hs_add_gene_symbol_from_ensembl_ids(colname_geneids_from="gene", colname_geneids_to="gene_name") %>% # map gene
    select(gene_name, everything(), -gene)
  df.es_mean
  
  df.es_mean.gene_filter <- df.es_mean %>% filter(gene_name %in% genes_select)
  df.es_mean.gene_filter <- df.es_mean.gene_filter %>% gather(key="annotation", value=es_weight, -gene_name)
  
  return(df.es_mean.gene_filter)
  # gene_name annotation es_weight
  # AGRP      ABC                0
  # POMC      ABC                0
  # AGRP      ACBG               0
  # POMC      ACBG               0
}

es_plot_gene_centric <- function(df, annotations_highlight, show_only_nonzero_es=F) {
  # this function supports multiple genes (genes_select as vector)
  # this function supports multiple annotations (annotations_highlight as vector)
  
  ### For DEV
  # df <- df.es
  # annotations_highlight <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12") # BMI_UKBB_Loh2018 FDR sign cell-types mousebrain
  # show_only_nonzero_es <- F
  
  if (show_only_nonzero_es) {
    df <- df %>% filter(es_weight > 0)  
  }
  
  ### "Ordering categories within ggplot2 facets" 
  ### REF: https://drsimonj.svbtle.com/ordering-categories-within-ggplot2-facets
  df <- df %>%
    # 2. Arrange by
    #   i.  facet group (gene_name)
    #   ii. bar height (es_weight)
    arrange(gene_name, es_weight) %>%
    # 3. Add order column of row numbers
    mutate(x_order = row_number())
  
  ### Plot
  p <- ggplot(df, aes(x=x_order, y=es_weight, label=annotation)) + 
    geom_point(size=2, color="gray") + 
    # geom_point(data=df %>% filter(annotation %in% annotations_highlight), color="black", size=3) +  # highlight black
    geom_point(data=df %>% filter(annotation %in% annotations_highlight), aes(color=es_weight), size=3) +  # highlight color scale
    scale_color_viridis_c(limits=c(0,1), direction=-1) + 
    geom_segment(aes(x=x_order, xend=x_order, y=0, yend=es_weight), color="grey") +
    geom_label_repel(data=df %>% filter(annotation %in% annotations_highlight)) +
    facet_wrap(~gene_name, 
               scales = "free" # x scale must be set 'free' because we use x='x_order', so every position is different
    ) +
    ### Axis
    labs(x="Cell-type", y="mean nES weight") + 
    ### Add categories to axis | replace the numeric values (x_order) on each x-axis with the appropriate label
    scale_x_continuous(
      breaks = df$x_order,
      labels = df$annotation,
      expand = c(0,0) # this is needed (for some reason) to make the x-axis labels properly aligned to the data
    ) +
    scale_y_continuous(limits=c(0,1)) + # fix the scale of the plot
    theme_classic() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    theme(axis.text.x=element_text(size = rel(0.5))) # smaller x-axis labels
  # theme(axis.text.x=element_blank()) # hide x-axis labels
  p
  return(p)
}


.es_plot_gene_centric_pathwork <- function(df, annotations_highlight, show_only_nonzero_es=F) {
  # TODO: make this function wrap plots using patchwork. 1: make a list of plots. 2: combine list of plots using patchwork::wrap_plots().
  stop("Not yet implemented")
  
  ### Plot [THIS CODE WORKS - save for now]
  p <- ggplot(df, aes(x=fct_reorder(annotation, es_weight), y=es_weight, label=annotation)) + 
    geom_point(size=3, color="gray") + 
    geom_point(data=df %>% filter(annotation %in% annotations_highlight), color="black") + 
    geom_segment(aes(x=annotation, xend=annotation, y=0, yend=es_weight), color="grey") +
    geom_label_repel(data=df %>% filter(annotation %in% annotations_highlight)) + # OK, but easier to use gghighlight
    # gghighlight(es_weight > 0, label_key=annotation) +  # OK, but a but too much labels when the x-axis is ordered by es_weight (the labels clutter)
    # gghighlight(annotation %in% annotations_highlight, label_key=NULL) + # WORKS - no label
    # gghighlight(annotation %in% annotations_highlight, label_key=annotation) + # WORKS - labels [*USED PRESENTATION*] [but does not support annotations_highlight as vector]
    # geom_label(data=df %>% filter(annotation %in% annotations_highlight)) + # give same results as gg_repel |  vjust="inward",hjust="inward"
    # gghighlight(annotation %in% annotations_highlight, label_key=annotation, label_params = list(hjust=0, vjust=0, direction="y")) + # cannot make 'label_params = list(hjust=0, vjust=0, direction="y")' to really do anything
    # gghighlight(annotation %in% annotations_highlight, label_key=annotation, label_params = list(nudge_x=1, nudge_y=1, position=NULL)) + # Error: Specify either `position` or `nudge_x`/`nudge_y`
    labs(x="Cell-type", y="mean nES weight") + 
    # scale_color_viridis("plasma", direction=-1) +
    # coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
    scale_y_continuous(limits=c(0,1)) # fix the scale of the plot
  # theme(plot.margin = unit(c(1,3,1,1), "cm")) # t, r, b, l 
  # theme(axis.text.x=element_blank()) # hide x-axis labels
  p
}
