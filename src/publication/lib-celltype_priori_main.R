############### SYNOPSIS ###################
# Helper functions for cell-type prioritization main figures


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
# =========================== READ h2 data and preprocess it ============ #
# ======================================================================= #

get_h2_data_with_appropriate_liability_scale <- function() {
  print("This function ensures that liability scale h2 estimates REPLACES observed scale h2 estimate for appropriate gwas'es (_liability_scale).")
  ### Load h2 estimates
  # All traits have observed scale h2 calculated
  # Some traits have liability scaled h2 calculated
  file.h2 <- here("results/h2_trait.multi_gwas.csv")
  df.h2 <- read_csv(file.h2)
  
  liability_pattern <- "_liability_scale"
  gwas.with_liability_scale <- df.h2 %>% filter(grepl(liability_pattern, gwas)) %>% pull(gwas)
  gwas.to_exclude <- str_replace_all(gwas.with_liability_scale, liability_pattern, "") # get 'basename' of liability gwas
  print("Will replace observed scale h2 with liability scale h2 for the following gwas'es:")
  print(gwas.to_exclude)
  ### Exclude gwas for which we have liability scale:
  df.h2 <- df.h2 %>% filter(!gwas %in% gwas.to_exclude)
  ### Rename liability gwas to 'basename'
  df.h2 <- df.h2 %>% mutate(gwas=str_replace_all(gwas, liability_pattern, ""))
  return(df.h2)
}



# ======================================================================= #
# ============ Calculate tau normalized (LOAD and ADD GWAS h2) ========== #
# ======================================================================= #

# M_NUMBER_SNPS <- 5961159 # EUR 1000G, MAF>=5%
### reference SNPs for European LD score estimation were the set of 9,997,231 SNPs with a minor allele count >= 5 from 489 unrelated European individuals in Phase 3 of 1000 Genomes.  ===> This is the number of SNPs in the annot files and in the .bim/.bed files. These are the SNPs used to calculate the ldscores (only some of the SNPs are written because of the --print-snps argument).
### Heritability was partitioned for the set of 5,961,159 SNPs with MAF >= 0.05 ===> This is the number of SNPs in the LDSC baseline model (number of SNPs that have written ldscores for)
### Regression coefficient estimation was performed with 1,217,312 HapMap3 SNPs (SNPs in HapMap3 are used because they are generally well-imputed). ===> regression weights
### SEE ALSO EVERNOTE "SOFTWARE | LDSC - understanding | SNPs in LDSC (weights, ldscores, summarystats...)"

### Calculate tau norm
### df.ldsc_cts <- df.ldsc_cts %>% left_join(df.h2, by="gwas") # join
### normalized tau = tau/(h2g/M), where h2g/M is the"mean per-SNP heritability"
### [LDSD-SEG]:    The normalized tau can be interpreted as the proportion by which the per-SNP heritability of an average SNP would increase if tau_k were added to it.
### [Benchmarker]: The normalized tau can be interpreted as the change in heritability associated with the value of the annotation increasing from 0 to 1.
# df.ldsc_cts <- df.ldsc_cts %>% mutate(tau_norm = estimate/(h2/M_NUMBER_SNPS))



# ======================================================================= #
# ============================ *UTIL FUNCTIONS* ========================= #
# ======================================================================= #

get_gwas_space_block <- function(gwas_id_vector) {
  ### DESCRIPTION: map vector from gwas_id to the spacing block
  
  ### DEV
  # gwas_id_vector <- "BMI_UKBB_Loh2018"
  
  block1 <- c("BMI_UKBB_Loh2018")
  block2 <- c("EA3_Lee2018",
              "INTELLIGENCE_Savage2018",
              "SCZ_Pardinas2018",
              "INSOMNIA_Jansen2018")
  block3 <- c("MS_Patsopoulos2011",
              "RA_Okada2014")
  block4 <- c("LIPIDS_LDL_Teslovich2010",
              "WHRadjBMI_UKBB_Loh2018",
              "HEIGHT_UKBB_Loh2018")
  out <- case_when(gwas_id_vector %in% block1 ~ 1,
                   gwas_id_vector %in% block2 ~ 2,
                   gwas_id_vector %in% block3 ~ 3,
                   gwas_id_vector %in% block4 ~ 4,
                   TRUE ~ NA_real_) # default value
  return(out)
}


get_gwas_filter <- function() {
  ### DESCRIPTION: returns ordered gwas traits to use for plotting
  
  ### Potential GWAS to include
  # "AN_PGC_Duncan2017"
  # "LUPUS_2015"
  filter.gwas <- c("BMI_UKBB_Loh2018",
                   #"NULL1",
                   "EA3_Lee2018",
                   "INTELLIGENCE_Savage2018",
                   "SCZ_Pardinas2018",
                   "INSOMNIA_Jansen2018",
                   #"NULL2",
                   "MS_Patsopoulos2011",
                   "RA_Okada2014",
                   #"NULL3",
                   "LIPIDS_LDL_Teslovich2010", # "LIPIDS_LDL_Willer2013",
                   "WHRadjBMI_UKBB_Loh2018",
                   "HEIGHT_UKBB_Loh2018"
  )
  
  # ### GWAS that should be on liability threshold
  # liability_gwas <- c("SCZ_Pardinas2018",
  #                     "INSOMNIA_Jansen2018",
  #                     "MS_Patsopoulos2011",
  #                     "RA_Okada2014"
  # )
  # ### append liability name
  # if (with_liability_names) {
  #   filter.gwas <- if_else(filter.gwas %in% liability_gwas, paste0(filter.gwas, "_liability_scale"), filter.gwas)
  # }
  return(filter.gwas)
}


# ======================================================================= #
# ============================ MULTI-GWAS PLOT =============================== #
# ======================================================================= #

set_multi_gwas_heatmap_plot <- function(df.ldsc_cts) {
  ### Create 'plotting' data frame
  filter.gwas <- get_gwas_filter() # get filters
  df.plot.multi_gwas <- df.ldsc_cts %>% filter(gwas %in% filter.gwas)
  ### Add "spacing block" column
  df.plot.multi_gwas <- df.plot.multi_gwas %>% mutate(space_block=get_gwas_space_block(gwas))
  
  ### Order GWAS
  df.plot.multi_gwas <- df.plot.multi_gwas %>% mutate(gwas=factor(gwas, levels=filter.gwas)) # Order GWAS
  df.plot.multi_gwas
  
  ### Rename GWAS [do this after filtering, so ensure uniqueness]
  ### IMPORANTLY recode_factor() allows us to rename the factor values but KEEP the ORDER defined by filter.gwas
  tmp_gwas_vector <- filter.gwas # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
  newnames <- utils.rename_gwas(tmp_gwas_vector, style="fullname") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
  rename_vector <- newnames; names(rename_vector) <- tmp_gwas_vector
  df.plot.multi_gwas <- df.plot.multi_gwas %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html
  ### fct_recode(): named character vectors where the name gives the new level, and the value gives the old level. 
  ### [did not work 1 - fct_recode] | newnames <- utils.rename_gwas(df.plot.multi_gwas$gwas, style="fullname") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
  ### [did not work 1 - fct_recode] | rename_vector <- as.character(df.plot.multi_gwas$gwas); names(rename_vector) <- newnames
  ### [did not work 1 - fct_recode] | df.plot.multi_gwas <- df.plot.multi_gwas %>% mutate(gwas_fmt = fct_recode(rename_vector))
  ### [did not work 2 - fct_relabel] | fct_relabel(unique(df.plot.multi_gwas$gwas), utils.rename_gwas, style="fullname") # does not help us becuase it works on the factor itself (unique values) and hence cannot be used with mutate()
  
  
  ### Plot
  p.multi_gwas <- ggplot(df.plot.multi_gwas, aes(x=gwas_fmt, y=annotation, fill=-log10(p.value))) + 
    geom_tile() + 
    geom_text(data=df.plot.multi_gwas %>% filter(fdr_significant), label="*", color="black", hjust=0.5, vjust=0.75) + # add asterisk if fdr significant
    # ^ hjust/vjust: https://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot
    # ^ hjust="center", vjust="middle"
    # scale_fill_viridis_c(option="magma", direction=-1) + 
    # scale_fill_distiller(palette="Greys", direction=1) + # "Greys" "Blues" "Greens" "BluGn" "Reds"
    # [GOOD] scale_fill_distiller(palette="Blues", direction=1, limits=c(0,5), oob=scales::squish) + # "Greys" "Blues" "Greens" "BluGn" "Reds"
    # scale_fill_distiller(palette="Greys", direction=1, limits=c(0,1.8), na.value = "white") + # tau norm plot
    colorspace::scale_fill_continuous_sequential(palette="Blues 2", rev=TRUE, limits=c(0,5), oob=scales::squish, na.value = "white") + # "Blues 2","Blues 3","Purples 3"
    # ^ oob: Function that handles limits outside of the scale limits (out of bounds). scales::squish "squishes" values into the range specified by limits
    labs(fill=expression(-log[10](P[S-LDSC]))) +
    theme_minimal() + # set this first
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank()) + # remove grid
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    theme(legend.position="bottom") +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y = element_blank(),
          axis.title = element_blank())
  
  ### Add spacing to heatmap
  # REF: https://stackoverflow.com/questions/40156061/white-space-between-tiles-in-heatplot-ggplot
  ### Approaches
  # 1: use facet_wrap ---> BEST
  # 2: add "NULL" values (+ order factor levels)
  # 3: patchwork / cowplot::plot_grid
  # 4: add separator line : http://www.roymfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r/ ---> only heatmap2
  p.multi_gwas <- p.multi_gwas + 
    facet_grid(~space_block, scales='free_x', space="free_x") + 
    theme(strip.background = element_blank(), strip.text.x = element_blank()) # remov facet labels completely
  p.multi_gwas

  ### Make global
  # "<<-" to set global variables | REF: https://stackoverflow.com/a/1236721
  p.multi_gwas <<- p.multi_gwas
  print("setting global vars: p.multi_gwas")
}


# ======================================================================= #
# ================================ h2 barplot =========================== #
# ======================================================================= #



set_h2_barplot <- function() {
  df.h2 <- get_h2_data_with_appropriate_liability_scale()
  
  ### Create 'plotting' data frame
  filter.gwas <- get_gwas_filter() # get filters
  df.plot.h2 <- df.h2 %>% filter(gwas %in% filter.gwas)
  ### Add "spacing block" column
  df.plot.h2 <- df.plot.h2 %>% mutate(space_block=get_gwas_space_block(gwas))
  ### Order GWAS
  df.plot.h2 <- df.plot.h2 %>% mutate(gwas=factor(gwas, levels=filter.gwas)) # Order GWAS
  ### Rename GWAS [not needed because names are not shown]
  # newnames <- utils.rename_gwas(filter.gwas, style="fullname") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
  # rename_vector <- newnames; names(rename_vector) <- filter.gwas
  # df.plot.multi_gwas <- df.plot.multi_gwas %>% mutate(gwas_fmt = recode_factor(gwas, !!!rename_vector)) # Use a named character vector to recode factors with unquote splicing. | REF: https://dplyr.tidyverse.org/reference/recode.html
  
  ### Plot
  p.h2 <- ggplot(df.plot.h2, aes(x=gwas, y=h2)) + 
    geom_col(fill="gray") +
    labs(y=expression(h[S-LDSC]^{2})) +
    theme(axis.text.x=element_text(angle=45, hjust=1), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_blank())
  p.h2
  
  ### Add spacing using facet_wrap
  p.h2 <- p.h2 + 
    facet_grid(~space_block, scales='free_x', space="free_x") + 
    theme(strip.background = element_blank(), strip.text.x = element_blank()) # remov facet labels completely
  
  ### Remove x-text
  p.h2_no_x <- p.h2 + theme(axis.text.x=element_blank())

  ### Make blank barplot
  p.h2.blank <- ggplot(df.plot.h2, aes(x=gwas, y=h2)) + 
    geom_blank() +
    theme(panel.grid = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.background = element_blank())
  
  ### Make global
  # "<<-" to set global variables | REF: https://stackoverflow.com/a/1236721
  p.h2 <<- p.h2
  p.h2_no_x <<- p.h2_no_x
  p.h2.blank <<- p.h2.blank
  print("setting global vars: p.h2, p.h2_no_x, p.h2.blank")
}






# ======================================================================= #
# ============================ MOUSEBRAIN UTILS ========================= #
# ======================================================================= #


# ============================ TAXONOMY TEXT ========================= #

get_celltype_taxonomy_text_position.mb <- function(df.metadata, df.plot) {
  ### INPUT: df.metadata with the columns
  # - TaxonomyRank4_reduced1
  # - tax_order_idx_mb_fig1c
  ### df.plot is only parsed as an extra safety check that annotations will be ordered correctly on the plot
  
  order.as_is <- df.plot$annotation %>% levels()
  order.as_should_be <- with(df.plot, annotation[order(tax_order_idx_mb_fig1c, as.character(annotation))])
  
  if (!all(order.as_is==order.as_should_be)) {
    print("Annotations should be ordered by taxonomy (from MB Fig1c).")
    print("This is how the df.plot should be ordered (secondary order by their annotation name not needed):")
    print("df.plot <- df.plot %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tax_order_idx_mb_fig1c, as.character(annotation))])))")
    stop("Error: df.plot annotation column is not ordered correctly by TaxonomyRank4_reduced1")
  }
  
  ### Create data frame with taxonomy text data
  df.tax_text_position <- df.metadata %>% 
    group_by(TaxonomyRank4_reduced1) %>% 
    summarize(n_annotations_in_tax=n()) %>% # count how many annotations in each tax
    left_join(df.metadata %>% select(TaxonomyRank4_reduced1, tax_order_idx_mb_fig1c) %>% distinct(), by="TaxonomyRank4_reduced1") %>% # add 'tax_order_idx_mb_fig1c' to be able to sort factor correctly.
    # ^ OBS: we use 'distinct()' because df.metadata contains (intentional) duplicated rows for some combinations of {TaxonomyRank4_reduced1, tax_order_idx_mb_fig1c}.
    mutate(TaxonomyRank4_reduced1 = factor(TaxonomyRank4_reduced1, levels=TaxonomyRank4_reduced1[order(tax_order_idx_mb_fig1c)])) # IMPORTANT: order factor by the SAME ORDER as 'annotation' column is ordered by in df.plot
  
  ### Add information of the text's position in the plot
  df.tax_text_position <- df.tax_text_position %>% 
    arrange(TaxonomyRank4_reduced1) %>%
    mutate(pos_mid=cumsum(n_annotations_in_tax)-n_annotations_in_tax/2,
           pos_start=cumsum(n_annotations_in_tax)-n_annotations_in_tax,
           pos_end=cumsum(n_annotations_in_tax),
           idx=1:n()) # idx is used for identifying every n'th tax
  df.tax_text_position <- df.tax_text_position %>% mutate(flag_draw_rect=if_else(idx %% 2 == 0, TRUE, FALSE)) 
  
  ### Remove some text because they contain too few cell-types
  df.tax_text_position %>% arrange(n_annotations_in_tax) # ---> potentially filter : n_annotations_in_tax >= 5
  df.tax_text_position <- df.tax_text_position %>% mutate(TaxonomyRank4_reduced1 = case_when(
    TaxonomyRank4_reduced1 == "Other CNS neurons" ~ "",
    # TaxonomyRank4_reduced1 == "Other glia" ~ "",
    TRUE ~ as.character(TaxonomyRank4_reduced1))
  )
  # TaxonomyRank4_reduced1     n_annotations_in_tax tax_order_idx_mb_fig1c pos_mid pos_start pos_end   idx flag_draw_rect
  # 1 Immune/Microglia                              5                      1     2.5         0       5     1 FALSE         
  # 2 Vascular                                     10                      2    10           5      15     2 TRUE 
  return(df.tax_text_position)
}

# ============================ TAXONOMY TEXT ========================= #
get_celltype_priori_base_tax_plot.mb <- function(df.plot, df.tax_text_position) {
  ### INPUT: df.plot should contain the following columns:
  # annotation
  # p.value.mlog10
  fdr_threshold <- 0.05/nrow(df.plot)
  p.main <- ggplot() +
    ### add ggplot 'baselayer'. This makes our 'canvas' and needs to go first (for some reason I don't fully understand...)
    geom_blank(data=df.plot, aes(x=annotation, y=p.value.mlog10)) +
    ### add tax info and gray rects (add gray rect before the data points, so avoid them 'covering' the points)
    geom_text(data=df.tax_text_position, aes(x=pos_mid, y=-0.5, label=TaxonomyRank4_reduced1), hjust="right", size=rel(3)) +
    geom_rect(data=df.tax_text_position %>% filter(flag_draw_rect), aes(xmin=pos_start, xmax=pos_end, ymin=-3, ymax=Inf), color="gray", alpha=0.1) +
    ### cell-types
    geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
    ### axes
    labs(x="", y=expression(-log[10](P[S-LDSC]))) +
    # coord
    coord_flip(ylim = c( 0, max(df.plot$p.value.mlog10) ), # This focuses the y-axis on the range of interest
               clip = 'off') +   # This keeps the labels from disappearing 
    # ^ clip = 'off': it allows drawing of data points anywhere on the plot, including in the plot margins. If limits are set via xlim and ylim and some data points fall outside those limits, then those data points may show up in places such as the axes, the legend, the plot title, or the plot margins.
    # ^ clip = 'off': disable cliping so df.tax_text_position text annotation can be outside of plot | REF https://stackoverflow.com/a/51312611/6639640
    ### theme
    theme_classic() + 
    theme(axis.text.y=element_text(size=rel(0.2)),
          axis.ticks.y=element_blank())
  p.main <- p.main + theme(plot.margin = unit(c(1,1,1,10), "cm")) # (t, r, b, l) widen left margin
  return(p.main)
}







