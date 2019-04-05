
# ======================================================================= #
# ===================== REQUIREMENTS FOR THIS LIB =================== #
# ======================================================================= #

# df.ldsc_cts


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
  newnames <- utils.rename_gwas(filter.gwas, style="fullname") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
  rename_vector <- newnames; names(rename_vector) <- filter.gwas
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
    labs(fill=expression(-log[10](P))) +
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
    labs(y=expression(h[SLDSC]^{2})) +
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


