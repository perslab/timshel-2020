

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
  # annotations <- c("TEGLU23","DEINH3")
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
# ================================ ES TABLES ================================ #
# ==================================================================================== #

get_annotation_es.top_n <- function(sem_obj, annotations, n_top_genes, es_metric) {
  # returns data frame with top n ES genes for a given metric
  # output data frame will have the column "es_genes_top_fmt"
  ### DEV
  # annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
  # n_top_genes <- 10
  # es_metric <- "es_mu"
  
  .check_valid_es_metrics(es_metric)
  stopifnot(length(es_metric)==1) # only one es metric allowed
  
  
  list.res <- list()
  for (annotation in annotations) {
    print(sprintf("Fetching data for annotation = %s", annotation))
    df.es <- get_es.annotation_centric(sem_obj, annotation=annotation) %>% # TODO: suppress print messages
      filter(es_metric == (!!es_metric) ) # !!sym(es_metric) does not work here [Error: object 'es_mu' not found]. REF: https://github.com/tidyverse/dplyr/issues/3139
    # gene            gene_name gene_idx_fixed es_metric es_weight gene_idx_sorted_within_metric
    # <chr>           <chr>              <int> <chr>         <dbl>                         <int>
    # 1 ENSG00000196344 ADH7                1319 es_mu         0.999                         15071
    # 2 ENSG00000171540 OTP                 1322 es_mu         0.998                         15070
    # 3 ENSG00000159723 AGRP                1895 es_mu         0.997                         15069
    df.es <- df.es %>% arrange(desc(es_weight)) %>% slice(1:n_top_genes)
    df.es <- df.es %>% select(gene, gene_name, es_weight)
    list.res[[annotation]] <- df.es
  }
  df.res <- bind_rows(list.res, .id="annotation")
  df.res <- df.res %>% 
    group_by(annotation) %>%
    arrange(desc(es_weight)) %>% 
    mutate(es_genes_top_fmt = paste0(gene_name, collapse=", ")) %>% # genes are ordered by es_weight so they will appear in the correct order
    ungroup()
  # if (!return_full_results) {
  #   df.res <- df.res %>% select(annotation, es_genes_top_fmt) %>% distinct()
  # }
  return(df.res)
}

get_annotation_es.table <- function(sem_obj, annotations, es_metric) {
  # returns data frame: genes x annotation
  ### DEV
  # annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
  # n_top_genes <- 10
  # es_metric <- "es_mu"
  
  .check_valid_es_metrics(es_metric)
  stopifnot(length(es_metric)==1) # only one es metric allowed
  
  
  list.res <- list()
  for (annotation in annotations) {
    print(sprintf("Fetching data for annotation = %s", annotation))
    df.es <- get_es.annotation_centric(sem_obj, annotation=annotation) %>% # TODO: suppress print messages
      filter(es_metric == (!!es_metric) ) # !!sym(es_metric) does not work here [Error: object 'es_mu' not found]. REF: https://github.com/tidyverse/dplyr/issues/3139
    # gene            gene_name gene_idx_fixed es_metric es_weight gene_idx_sorted_within_metric
    # <chr>           <chr>              <int> <chr>         <dbl>                         <int>
    # 1 ENSG00000196344 ADH7                1319 es_mu         0.999                         15071
    # 2 ENSG00000171540 OTP                 1322 es_mu         0.998                         15070
    # 3 ENSG00000159723 AGRP                1895 es_mu         0.997                         15069
    df.es <- df.es %>% select(gene, gene_name, es_weight)
    list.res[[annotation]] <- df.es
  }
  df.res <- bind_rows(list.res, .id="annotation")
  # annotation gene            gene_name es_weight
  # 1 TEGLU23    ENSG00000141668 CBLN2         0    
  # 2 TEGLU23    ENSG00000204624 DISP3         0.228
  # 3 TEGLU23    ENSG00000187848 P2RX2         0    
  df.res.spread <- df.res %>% spread(key="annotation", value="es_weight")
  # gene            gene_name DEGLU4 DEGLU5  DEINH3 MEGLU1 MEGLU10 MEGLU11 MEINH2  TEGLU17 TEGLU23 TEGLU4 TEINH12
  # 1 ENSG00000000003 TSPAN6     0     0      0        0       0       0     0      0        0       0       0     
  # 2 ENSG00000000005 TNMD       0.501 0.390  0.127    0.881   0       0     0.495  0        0       0.857   0.105 
  # 3 ENSG00000000457 SCYL3      0.102 0      0.00221  0       0       0     0.0304 0.168    0.254   0.0362  0     
  return(df.res.spread)
}


# ==================================================================================== #
# ================================ Annotation centric ================================ #
# ==================================================================================== #
### Annotation centric: 
# Plot ES of all genes for **one annotation** (annotation). 
# You can highlight genes in the plot function.

get_es.annotation_centric <- function(sem_obj, annotation) {
  ### OBS: this function only supports extracting one annotation.
  
  if (is.null(sem_obj[["group_by_annotation.sem"]][[annotation]])) {
    stop(sprintf("Annotation = '%s' not found in data", annotation))
  }
  ### Get data
  df.es <- sem_obj[["group_by_annotation.sem"]][[annotation]] %>% 
    mutate(gene=sem_obj$genes) %>% 
    hs_add_gene_symbol_from_ensembl_ids(colname_geneids_from="gene", colname_geneids_to="gene_name") # map gene
  ### Add addition data
  df.es <- df.es %>% 
    mutate(
      expr_mean = sem_obj$data$mean[[annotation]],
      es_mu = sem_obj$sem_meta$mean[[annotation]] 
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
  #   ENSG00000141668 CBLN2                  1 tstat    -0.885                          1040
  #   ENSG00000204624 DISP3                  2 tstat    -0.653                          7517
  #   ENSG00000187848 P2RX2                  3 tstat    -1.33                            211
  #   ENSG00000171522 PTGER4                 4 tstat    -0.722                          4716
}


plot_es.annotation_centric <- function(df.es, genes_highlight) {
  # df.es <- df.es %>% mutate(gene_label = paste0(gene_name, "(", round(es_weight,2), ")") ) # looks too complication
  # TODO: try to add geom_density2d() to the plot
  
  ### rename
  df.es <- df.es %>% mutate(es_metric = case_when(
    es_metric == "es_mu" ~ expression(ES[mu]),
    es_metric == "expr_mean" ~ "Mean norm. expr.",
    es_metric == "ges" ~ "GES",
    es_metric == "specificity" ~ "EP",
    es_metric == "si" ~ "NSI",
    es_metric == "tstat" ~ "DET",
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
# ================================ UTILS ================================ #
# ============================================================================== #

.get_genes_select_idx <- function(sem_obj, genes_select) {
  ### Function to get indices of genes in genes_select. Genes not in sem_obj will be dropped.
  ### Input: genes_select is a character vector of gene names
  ### Return: a named character vector. values are indecies of matching genes, names are gene names.
  
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
  names(genes_match.idx) <- genes_match.names # set names
  return(genes_match.idx)
}

.check_valid_es_metrics <- function(es_metric) {
  es_metrics_allowed <- c("expr_mean", "es_mu", "tstat", "ges", "si", "specificity")
  if (!any(es_metric %in% es_metrics_allowed)) {
    stop(sprintf("Argument for es_metric=%s did not match es_metrics_allowed: [%s]", es_metric, paste(es_metrics_allowed, collapse=", ")))
  }
}

# ============================================================================== #
# ================================ Gene centric ================================ #
# ============================================================================== #
### Gene centric: Plot ES of all annotations for **one gene** (genes_select supports multiple genes as subplots). 
# You can highlight annotations in the plot function.

get_es.gene_centric.all_es_metrics <- function(sem_obj, genes_select) {
  ### returns all ES metrics (incl. es_mu and mean expr)
  # this function supports multiple genes (genes_select as vector)
  
  ### DEVELOPMENT
  # sem_obj
  # genes_select <- c("AGRP", "NPY")
  
  genes_match.idx <- .get_genes_select_idx(sem_obj, genes_select)

  list.es <- list()
  for (annotation in sem_obj$annotations) {
    df.es <- sem_obj[["group_by_annotation.sem"]][[annotation]] %>% slice(genes_match.idx) # df | ges, si, tstat, specificity
    expr_mean <- sem_obj$data$mean %>% slice(genes_match.idx) %>% pull(!!sym(annotation)) # vector | mean expression
    es_mean <- sem_obj$sem_meta$mean %>% slice(genes_match.idx) %>% pull(!!sym(annotation)) # vector | es_mu
    list.es[[annotation]] <- df.es %>% mutate(
      expr_mean=expr_mean,
      es_mean=es_mean,
      gene_name=names(genes_match.idx))
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



get_es.gene_centric.single_es_metric <- function(sem_obj, genes_select, es_metric) {
  # this function supports multiple genes (genes_select as vector)
  
  .check_valid_es_metrics(es_metric)
  stopifnot(length(es_metric)==1) # only one es metric allowed

  
  ### Get matching genes
  genes_match.idx <- .get_genes_select_idx(sem_obj, genes_select)
  
  
  ### Get data
  if (es_metric == "es_mu") {
    df.es <- es_mean <- sem_obj$sem_meta$mean
  } else if (es_metric == "expr_mean") {
    df.es <- sem_obj$data$mean
  } else { # 'real' es metric
    df.es <- sem_obj$sem[[es_metric]]
  }
  
  df.es_mean.gene_filter <- df.es %>% 
    slice(genes_match.idx) %>% # select genes
    mutate(gene_name=names(genes_match.idx)) # add gene names

  ### OLD method | WORKS [SIMPLE] [but OK delete]
  # df.es_mean <- sem_obj$sem_meta$mean %>% 
  #   mutate(gene=sem_obj$genes) %>% 
  #   hs_add_gene_symbol_from_ensembl_ids(colname_geneids_from="gene", colname_geneids_to="gene_name") %>% # map gene
  #   select(gene_name, everything(), -gene)
  # df.es_mean.gene_filter <- df.es_mean %>% filter(gene_name %in% genes_select)
  
  ### Transform to long format
  df.es_mean.gene_filter <- df.es_mean.gene_filter %>% gather(key="annotation", value=es_weight, -gene_name)
  
  ### order genes by their order in genes_select [this makes it easier to get a nicely ordered facet plot in plot_es.gene_centric.single_es_metric()]
  df.es_mean.gene_filter <- df.es_mean.gene_filter %>% mutate(gene_name=factor(gene_name, levels=genes_select))
  
  return(df.es_mean.gene_filter)
  # gene_name annotation es_weight
  # AGRP      ABC                0
  # POMC      ABC                0
  # AGRP      ACBG               0
  # POMC      ACBG               0
}

plot_es.gene_centric.single_es_metric <- function(df, 
                                                  n_top_annotations_highlight=NULL, # integer specifying how many annotations to highlight
                                                  annotations_highlight=NULL, # character vector of annotations to highlight. If NULL, 
                                                  # if annotations_highlight is a named vector, annotations will be renamed
                                                  # using names(annotations_highlight) as new_names and the annotations_highlight to match the original names.
                                                  annotations_colormap=NULL, # named vector: x=colors, names(x)=annotations
                                                  show_only_nonzero_es=F, 
                                                  scale.zero_to_one=T, # y-scale (0-1 for es_mu)
                                                  scale.log=F, # y-scale log-transform
                                                  size.highlight.text=rel(3.5),
                                                  size.highlight.points=rel(2)
                                                  ) {
  # this function supports multiple genes (genes_select as vector)
  # this function supports multiple annotations (annotations_highlight as vector)
  ### INPUT df
  # gene_name annotation es_weight
  # 1 POMC      ABC           -2.03 
  # 2 AGRP      ABC           -0.616
  # 3 LEPR      ABC            1.71 
  
  ### For DEV
  # df <- df.es.gene
  # annotations_highlight <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12") # BMI_UKBB_Loh2018 FDR sign cell-types mousebrain
  # annotations_colormap <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")
  # n_top_annotations_highlight <- 3
  # show_only_nonzero_es <- F
  # scale.zero_to_one <- T
  # scale.log <- F

  ### Error checks
  if (!is.null(n_top_annotations_highlight) & !is.null(annotations_highlight)) {
    stop("You must decide between n_top_annotations_highlight and annotations_highlight.")
  }
  if (is.null(annotations_highlight) & !is.null(annotations_colormap)) {
    stop("You must supply annotations_highlight when using annotations_colormap.")
  }
  if (!is.null(annotations_highlight) & !is.null(annotations_colormap)) {
    if (length(annotations_highlight) != length(annotations_colormap)) {
      stop("annotations_highlight and annotations_colormap must have same length.")
      # This is to avoid downstream errors, but might not strictly be needed in all cases
    }
  }
  if (!is.null(n_top_annotations_highlight) & !is.null(annotations_colormap)) {
    stop("You must decide between n_top_annotations_highlight and annotations_highlight: only supply one or the arguments.")
  }
  

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
  # gene_name annotation es_weight x_order
  # 1 AGRP      MOL2           -4.37       1
  # 2 AGRP      ACNT1          -3.75       2
  # 3 AGRP      MOL1           -3.66       3
  # 
  ### Set flag_highlight
  if (!is.null(annotations_highlight)) {
    if (!all(annotations_highlight %in% df$annotation)) {
      print("*WARNING*: not all annotations_highlight were found in the input data frame.")
    }
    df <- df %>% mutate(flag_highlight = if_else(annotation %in% annotations_highlight, TRUE, FALSE))
  } else {
    df <- df %>% 
      group_by(gene_name) %>%
      mutate(flag_highlight = if_else(rank(-es_weight) <= n_top_annotations_highlight, TRUE, FALSE)) %>%
      ungroup()
  }
  
  if (rlang::is_named(annotations_highlight)) { # ALT: if (!is.null(names(annotations_highlight))) {}
    print("Received named annotations_highlight. Will rename annotations (incl. annotations_colormap if given)")
    rename_vector <- names(annotations_highlight)
    names(rename_vector) <- annotations_highlight
    df <- df %>% mutate(annotation = recode(annotation, !!!rename_vector)) # x=new_names; names(x)=old_names
    if (!is.null(annotations_colormap)) {
      names(annotations_colormap) <- recode(names(annotations_colormap), !!!rename_vector)
    }
  }
  
  
  ### Init plot
  p <- ggplot(df, aes(x=x_order, y=es_weight, label=annotation)) + 
    geom_segment(aes(x=x_order, xend=x_order, y=0, yend=es_weight), color="grey", alpha=0.3) +
    geom_point(size=size.highlight.points, color="gray")
  ### colors
  if (!is.null(annotations_colormap)) {
    p <- p + 
      geom_point(data=df %>% filter(flag_highlight), aes(color=annotation), size=size.highlight.points) +
      geom_text_repel(data=df %>% filter(flag_highlight), aes(color=annotation), size=size.highlight.text, segment.alpha=0.5)
  } else {
    p <- p +
      geom_point(data=df %>% filter(flag_highlight), color="dodgerblue", size=size.highlight.points) +  # highlight color scale
      geom_text_repel(data=df %>% filter(flag_highlight), color="dodgerblue", size=size.highlight.text, segment.alpha=0.5)
  }
  ### Theme
  p <- p +  
    theme_classic() + # 'base theme'
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    theme(axis.text.x=element_text(size = rel(0.15))) # smaller x-axis labels
  # theme(axis.text.x=element_blank()) # hide x-axis labels
  ### facet wrap
  if (n_distinct(df$gene_name) > 1) {
  ### Alternative to stringr::str_to_title
  # capitalize <- function(string) {
  #   substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  #   return(string)
  # }
  p <- p +  
    facet_wrap(~gene_name, 
               scales = "free", # x scale must be set 'free' because we use x='x_order', so every position is different
               labeller = as_labeller(stringr::str_to_title) # capitalize first letter
    ) + 
    theme(strip.background=element_blank())
  }
  ### Axis
  p <- p +  labs(x="Cell-type", y="ES weight")
  ### Add categories to axis | replace the numeric values (x_order) on each x-axis with the appropriate label
  p <- p +  scale_x_continuous(
    breaks = df$x_order,
    labels = df$annotation,
    expand = c(0,0) # this is needed (for some reason) to make the x-axis labels properly aligned to the data
  )
  ### Scale and colors
  if (scale.zero_to_one) {
    p <- p + scale_y_continuous(limits=c(0,1))
  }
  # if (scale.zero_to_one & is.null(annotations_colormap)) {
  #   p <- p + scale_color_viridis_c("ES weight", limits=c(0,1), direction=-1)
  # }
  if (!is.null(annotations_colormap)) {
    p <- p + scale_color_manual(values=annotations_colormap)
    p <- p + guides(color=F) # hide legend
  }
  if (scale.log) {
    # p <- p + scale_color_viridis_c("ES weight", direction=-1)
    p <- p + scale_y_log10()
  }
  return(p)
}



.es_plot_gene_centric_pathwork <- function(df, annotations_highlight, show_only_nonzero_es=F) {
  # TODO: make this function wrap plots using patchwork. 
  # 1: make a list of plots. 
  # 2: combine list of plots using patchwork::wrap_plots().
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
