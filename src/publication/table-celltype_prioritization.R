############### SYNOPSIS ###################
### AIM: Write cell-type prioritization result tables

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ================================= RUN ================================= #
# ======================================================================= #


dataset_prefixes <- c("tabula_muris", "mousebrain", "hypothalamus")
filter.gwas_bmi <- "BMI_UKBB_Loh2018"
filter.gwas <- utils.get_gwas_ids_for_meta_analysis()

for (dataset_prefix in dataset_prefixes) {
  print(dataset_prefix)
  # file.data <- here("results", sprintf("prioritization_celltypes--%s.multi_gwas.csv.gz", dataset_prefix))
  file.data <- here("results/cellect_ldsc/prioritization.csv")
  df <- suppressMessages(read_csv(file.data))
  df <- format_cellect_ldsc_results(df)
  # gwas               estimate     std.error p.value annotation       
  # 1 AD_Jansen2019 0.00000000395 0.00000000218  0.0348 n03              
  # 2 AD_Jansen2019 0.00000000625 0.00000000354  0.0389 s15.Microglia 
  if (dataset_prefix == "hypothalamus") {
    df <- df %>% filter(specificity_id %in% get_scrna_seq_dataset_prefixes("hypo"))
  } else {
    df <- df %>% filter(specificity_id == dataset_prefix)
  }
  # =================== BMI GWAS LDSC CTS ================== #
  df.bmi <- df %>% filter(gwas %in% filter.gwas_bmi) %>% select(-gwas)
  ### Add fdr_significant flag for BMI
  # df.ldsc_cts.tmp.summary <- df.ldsc_cts %>% group_by(specificity_id) %>% summarise(n_obs=n())
  # df.ldsc_cts <- left_join(df.ldsc_cts, df.ldsc_cts.tmp.summary, by="specificity_id")
  df.bmi <- df.bmi %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))
  # df.ldsc_cts
  
  # df.bmi <- df.bmi %>% select(annotation, "Coeficient P-value"=p.value, "Coefficient"=estimate, "Coefficient std. error"=std.error)
  # =================== PREP NON-BMI GWAS LDSC CTS P-values ================== #
  df.other_gwas <- df %>% filter(gwas %in% filter.gwas) # filter GWAS 'trait representative'
  df.other_gwas <- df.other_gwas %>% filter(!gwas %in% filter.gwas_bmi) # remove BMI results
  ### Rename GWAS
  tmp_gwas_vector <- filter.gwas # convenience selector. If you want a specific order of the result, this vector should contain (unique) ordered values
  newnames <- utils.rename_gwas(tmp_gwas_vector, style="abrv_author_year") # "fullname_author_year","fullname","abrv_author_year","abrv_year","abrv" 
  rename_vector <- newnames; names(rename_vector) <- tmp_gwas_vector
  df.other_gwas <- df.other_gwas %>% mutate(gwas = recode_factor(gwas, !!!rename_vector))
  ### spread
  df.other_gwas.spread <- df.other_gwas %>% 
    select(gwas, p.value, annotation, specificity_id) %>%
    spread(key="gwas", value="p.value")
  
  if (dataset_prefix=="mousebrain") { # only do h2 estimates for mousebrain
    # =================== BMI GWAS h2 results ================== #
    # h2 results exists for all mousebrain annotations and some TM annotations. 
    # Because we do left_join, we get the approriate join.
    file.h2_annotation <- here("results/cellect_ldsc/heritability.csv")
    df.h2 <- suppressMessages(read_csv(file.h2_annotation))
    df.h2
    # gwas       specificity_id annotation      Prop._SNPs Prop._h2 Prop._h2_std_er… h2_enrichment h2_enrichment_se h2_enrichment_p…
    # RA_Okada2… mousebrain     TEINH12              0.171    0.241          0.0158           1.41           0.0926         5.13e- 5
    ### OLD
    # run_name annotation gwas  Prop._SNPs Prop._h2 Prop._h2_std_er… Enrichment Enrichment_std_… Enrichment_p Coefficient
    # 1 celltyp… DEGLU4     AD_J…     0.0985    0.151           0.0217       1.54            0.221     0.00889    -1.32e- 9
    # 2 celltyp… DEGLU5     AD_J…     0.117     0.344           0.156        2.93            1.33      0.120       4.99e- 9
    # df.h2 <- df.h2 %>% select(-run_name, starts_with("Coefficient"))
    rename_select_vector <- c(
      "annotation", 
      "gwas",
      "Annotation size"="Prop._SNPs",
      "h2g prop."="Prop._h2",
      "h2g enrichment"="h2_enrichment",
      "h2g enrichment P-value"="h2_enrichment_pvalue") # new_name = old_name
    df.h2 <- df.h2 %>% select(!!!rename_select_vector)
    df.h2 <- df.h2 %>% filter(gwas %in% filter.gwas_bmi) %>% select(-gwas)
    df.join <- df.bmi %>% left_join(df.h2, by="annotation") # *TMP disabled 2019-11*
    df.join <- df.join %>% mutate(annotation_fmt = annotation) # needed to match TM/Hypo format
  } else if (dataset_prefix == "hypothalamus") {
    df.join <- df.bmi %>% mutate(annotation_fmt = utils.rename_annotations.hypothalamus(annotation, specificity_id, check_all_matches=T))
    df.other_gwas.spread <- df.other_gwas.spread %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, specificity_id, check_all_matches=T))
  } else if (dataset_prefix == "tabula_muris") {
    df.join <- df.bmi %>% mutate(annotation_fmt = utils.rename_annotations.tabula_muris(annotation, style="tissue - celltype", check_all_matches=T))
  } else {
    stop("Internal error")
  }
  
    
  # =================== COMBINE ================== #
  df.join <- df.join %>% left_join(df.other_gwas.spread %>% select(-specificity_id), by=c("annotation_fmt"="annotation"))
  df.join <- df.join %>% arrange(p.value) # sort by BMI p-value
  df.join <- df.join %>% select("Annotation"=annotation_fmt, 
                                "Coefficient p-value"=p.value, 
                                "Coefficient"=estimate, 
                                "Coefficient std error"=std.error, 
                                "FDR < 0.05"=fdr_significant,
                                everything(),
                                -specificity_id) # rename and reorder cols
  file.out <- here("src/publication/tables", sprintf("table-celltype_ldsc_results.%s.csv", dataset_prefix))
  df.join %>% write_csv(file.out)
}





