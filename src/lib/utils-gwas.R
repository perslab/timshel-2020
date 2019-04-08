################ SYNOPSIS ###################
# Functions to rename GWAS trait names


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
# ===================== DATA GWAS_ID META-ANALYSIS ====================== #
# ======================================================================= #

# N_GWAS = 39
# GWAS represent represent the lastest GWAS for each trait.
# Traits are selected to be as independent/non-redundant as possible.
# LAST MODIFIED: 04.05.2019

.data.gwas_ids_meta_analysis <- "AD_Jansen2019
ADHD_PGC_Demontis2017
AN_PGC_Duncan2017
ASD_iPSYCH_PGC_Grove2018
BIP_PGC2018
BMI_UKBB_Loh2018
CAD_Schunkert2011
CARDIOVASCULAR_UKBB_Loh2018
CELIAC_Dubois2010
CROHNS_Jostins2012
DEPRESSION_Nagel2018
DIASTOLICadjMED_UKBB_Loh2018
EA3_Lee2018
FG_Female_Lagou2018
FI_Female_Lagou2018
HBA1C_MAGIC_Wheeler2017
HEIGHT_UKBB_Loh2018
IBD_Jostins2012
INSOMNIA_Jansen2018
INTELLIGENCE_Savage2018
LIPIDS_HDL_Teslovich2010
LIPIDS_LDL_Teslovich2010
LIPIDS_TG_Teslovich2010
LUPUS_Bentham2015
MDD_Howard2019
MS_Patsopoulos2011
NEUROTICISM_Nagel2018
PBC_Cordell2015
RA_Okada2014
RB_Linner_2019
SCZ_Pardinas2018
SWB_Okbay2016
SYSTOLICadjMED_UKBB_Loh2018
T1D_Bradfield2011
T2D_UKBB_Loh2018
UC_Jostins2012
WHR_Pulit2019
WHRadjBMI_UKBB_Loh2018
WORRY_Nagel2018"

# ======================================================================= #
# =========================== DATA RENAMING ============================= #
# ======================================================================= #
# TODO: move to .csv file in data/

.data.gwas_rename <- "trait_name_full	trait_name_abrv	author	year	gwas_id
Multiple Sclerosis	MS	Patsopoulos	2011	MS_Patsopoulos2011
Alzheimer's Disease	AD	Lambert	2013	AD_Lambert2013
Alzheimer's Disease	AD	Jansen	2019	AD_Jansen2019
Attention Deficit Hyperactivity Disorder	ADHD	Demontis	2017	ADHD_PGC_Demontis2017
Anorexia Nervosa 	AN	Duncan	2017	AN_PGC_Duncan2017
Autism Spectrum Disorder	ASD	Grove	2018	ASD_iPSYCH_PGC_Grove2018
Bipolar Disorder	BIP	PGC	2018	BIP_PGC2018
Major Depression Disorder	MDD	Howard	2019	MDD_Howard2019
Major Depression Disorder	MDD	Wray	2018	MDD_PGC_Wray2018
Schizophrenia	SCZ	Pardinas	2018	SCZ_Pardinas2018
Schizophrenia	SCZ	Ripke	2014	SCZ_EUR_Ripke2014
Educational Attainment	EA2	Okbay	2016	EA2_Okbay2016
Educational Attainment	EA3	Lee	2018	EA3_Lee2018
Insomnia	INS	Jansen	2018	INSOMNIA_Jansen2018
Intelligence	INT	Sniekers	2017	INTELLIGENCE_Sniekers2017
Intelligence	INT	Savage	2018	INTELLIGENCE_Savage2018
Risk Tolerance	RT	Linner	2019	RB_Linner_2019
Depressive Symptoms	DS	Okbay	2016	DS_Okbay2016
Neuroticism	NEU	Okbay	2016	NEUROTICISM_Okbay2016
Subjective Well Being	SWB	Okbay	2016	SWB_Okbay2016
Neuroticism	NEU	Nagel	2018	NEUROTICISM_Nagel2018
Depression	DEP	Nagel	2018	DEPRESSION_Nagel2018
Depressed Affect	DEP_AFFECT	Nagel	2018	DEPRESSED_AFFECT_Nagel2018
Worry	WORRY	Nagel	2018	WORRY_Nagel2018
BMI	BMI	Yengo	2018	BMI_UPDATE_Yengo2018
Height	HEIGHT	Yengo	2018	HEIGHT_Yengo2018
BMI	BMI	Locke	2015	BMI_Locke2015
Height	HEIGHT	Wood	2014	HEIGHT_Wood2014
Waist-hip Ratio	WHR	Shungin	2015	WHR_Shungin2015
Waist-hip Ratio (adj. BMI)	WHRadjBMI	Shungin	2015	WHRadjBMI_Shungin2015
Waist-hip Ratio	WHR	Pulit	2019	WHR_Pulit2019
Waist-hip Ratio (adj. BMI)	WHRadjBMI	Pulit	2019	WHRadjBMI_Pulit2019
BMI	BMI	Pulit	2019	BMI_Pulit2019
BMI (male)	BMI_M	Pulit	2019	BMI_male_Pulit2019
BMI (female)	BMI_F	Pulit	2019	BMI_female_Pulit2019
BMI	BMI	Loh	2018	BMI_UKBB_Loh2018
Waist-hip Ratio (adj. BMI)	WHRadjBMI	Loh	2018	WHRadjBMI_UKBB_Loh2018
Height	HEIGHT	Loh	2018	HEIGHT_UKBB_Loh2018
Type 2 Diabetes	T2D	Loh	2018	T2D_UKBB_Loh2018
Type 2 Diabetes (adj. BMI)	T2DadjBMI	Mahajan	2018	T2DadjBMI_DIAMANTE_Mahajan2018
Type 2 Diabetes	T2D	Mahajan	2018	T2D_UKBB_DIAMANTE_Mahajan2018
Type 2 Diabetes	T2D	Mahajan	2018	T2D_DIAMANTE_Mahajan2018
Type 2 Diabetes	T2D	Xue	2018	T2D_Xue2018
Fasting Insulin (male)	FI_M	Lagou	2018	FI_Male_Lagou2018
Fasting Glucose (male)	FG_M	Lagou	2018	FG_Male_Lagou2018
Fasting Insulin (female)	FI_F	Lagou	2018	FI_Female_Lagou2018
Fasting Glucose (female)	FG_F	Lagou	2018	FG_Female_Lagou2018
HbA1c	HBA1C	Wheeler	2017	HBA1C_MAGIC_Wheeler2017
High Density Lipoprotein	HDL	Willer	2013	LIPIDS_HDL_Willer2013
Low Density Lipoprotein	LDL	Willer	2013	LIPIDS_LDL_Willer2013
Triglycerides	TG	Willer	2013	LIPIDS_TG_Willer2013
Total Cholesterol	TC	Willer	2013	LIPIDS_TC_Willer2013
Type 1 Diabetes	T1D	Bradfield	2011	T1D_Bradfield2011
High Density Lipoprotein	HDL	Teslovich	2010	LIPIDS_HDL_Teslovich2010
Low Density Lipoprotein	LDL	Teslovich	2010	LIPIDS_LDL_Teslovich2010
Triglycerides	TG	Teslovich	2010	LIPIDS_TG_Teslovich2010
Celiac Disease	CELIAC	Dubois	2010	CELIAC_Dubois2010
Ulcerative Colitis	UC	Jostins	2012	UC_Jostins2012
Crohn's Disease	CROHNS	Jostins	2012	CROHNS_Jostins2012
Inflammatory Bowel Disease	IBD	Jostins	2012	IBD_Jostins2012
Primary Biliary Cirrhosis	PBC	Cordell	2015	PBC_Cordell2015
Systemic Lupus Erythematosus	SLE	Bentham	2015	LUPUS_Bentham2015
Rheumatoid Arthritis	RA	Okada	2014	RA_Okada2014
Coronary Artery Disease	CAD	Schunkert	2011	CAD_Schunkert2011
Diastolic Blood Pressure (adj. medicine)	BP_DIA	Loh	2018	DIASTOLICadjMED_UKBB_Loh2018
Systolic Blood Pressure (adj. medicine)	BP_SYS	Loh	2018	SYSTOLICadjMED_UKBB_Loh2018
Cardiovascular Diseases	CVD	Loh	2018	CARDIOVASCULAR_UKBB_Loh2018"


# ======================================================================= #
# ========================= FUNCTION: RENAME ============================= #
# ======================================================================= #

utils.get_gwas_ids_for_meta_analysis <- function() {
  gwas_ids.meta_analysis <- read_csv(.data.gwas_ids_meta_analysis, col_names=F) %>% pull(X1)
  print(sprintf("Returning n=%s gwas_ids", length(gwas_ids.meta_analysis)))
  return(gwas_ids.meta_analysis)
}
  

# ======================================================================= #
# ======================== FUNCTION: RENAME GWAS ========================= #
# ======================================================================= #

utils.rename_gwas <- function(gwas_ids, style, check_all_matches=F, return_as_df=F) {
  ### DESCRIPTION
  # This function is an alternative version of utils.rename_gwas.drop_no_matches()
  # gwas_id's not found in the data base, will not be renamed.
  # *Advantage of this function*: this function guarantees the same order and length of input-to-output.
  ### INPUT
  # gwas_ids:      a vector of gwas_ids's
  # style:        style.
  # check_all_matches:  if true, the function will raise an exception if the input gwas_ids are not all found in the database.
  ### OUTPUT
  # a vector of formatted gwas names. 
  # the length of the vector is equal to the input gwas_ids.
  
  ### DEV
  # style <- "fullname_author_year"
  # gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010", "asdasd") # must fail
  # gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010") # must pass
  # return_as_df <- FALSE
  # check_all_matches <- FALSE

  ### Check input
  allowed_styles <- c("fullname_author_year",
                      "fullname",
                      "abrv_author_year",
                      "abrv_year",
                      "abrv" # OBS: abrv may not be unique!
                      )
  if (!style %in% allowed_styles) {
    stop("Got wrong style. Allowed styles are: [%s]", paste(allowed_styles, collapse=", "))
  }
  if (style == "abrv") {
    print("Warning: trait name abreviations are not unique.")
  }
  ### Read data
  df.gwas_database <- read_tsv(.data.gwas_rename)
  if (any(duplicated(df.gwas_database$gwas_id))) {
    stop("internal renaming data frame contains duplicates. Fix the .csv file")
  }
  ### Check that all gwas_ids are present in "data base"
  bool.found <- gwas_ids %in% df.gwas_database$gwas_id
  if (check_all_matches && !all(bool.found)) {
    gwas_ids_not_found <- gwas_ids[!bool.found]
    stop(sprintf("Some gwas_ids not found: [%s]", paste(gwas_ids_not_found, collapse=",")))
  }
  ### Format output conditional on style
  df.gwas_database.select <- df.gwas_database %>% filter(gwas_id %in% gwas_ids)
  if (style=="fullname_author_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_full, " (", author, ", ", year, ")"))
  } else if (style=="fullname") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_full))
  } else if (style=="abrv_author_year") {
      df.gwas_out_fmt <- df.gwas_database.select %>% 
        mutate(fmt=paste0(trait_name_abrv, " (", author, ", ", year, ")"))
  } else if (style=="abrv_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_abrv, " (", year, ")"))
  } else if (style=="abrv") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_abrv))
  } else {
    stop("Internal function error")
  }
  ### Join renamed with input
  df.all <- tibble(gwas_id=gwas_ids) # this ensures that the return value has the same length as the input
  df.all <- df.all %>% left_join(df.gwas_out_fmt, by="gwas_id")
  ### Add fmt column for gwas_id non-matches
  ### df.all <- df.all %>% mutate(fmt=replace(fmt, is.na(fmt), gwas_id)) # GIVES WARNING ---> Warning message: In x[list] <- values :number of items to replace is not a multiple of replacement length
  df.all <- df.all %>% mutate(fmt=if_else(is.na(fmt), gwas_id, fmt)) 
  # ^ gwas_id non-matches have NA values in the 'fmt' column
  # ^ we replace na with the input gwas_id name
  if (return_as_df) {
    gwas_out_fmt <- df.all
  } else(
    gwas_out_fmt <- df.all %>% pull(fmt) # returns vector
  )
  return(gwas_out_fmt)
}

### EXAMPLE CALL
# style <- "fullname_author_year"
# gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010") # must pass
# gwas_ids.fmt <- utils.rename_gwas(gwas_ids, style, return_as_df=F)
# gwas_ids.fmt




# ======================================================================= #
# =============================== RENAME FUNCTION ================================= #
# ======================================================================= #
# THIS FUNCTION WORKS


utils.rename_gwas.drop_no_matches <- function(gwas_ids, style, check_all_matches=T, return_as_df=F) {
  ### INPUT
  # gwas_ids:           a vector of gwas_ids's
  # style:              style.
  # check_all_matches:  if true, the function will raise an exception if the input gwas_ids are not all found in the database.
  #                     if false, non-matches will be dropped silently.
  ### OUTPUT
  # a vector of formatted gwas names. The length is equal to the number of gwas_id that matched the database.
  
  ### DEV
  # style <- "fullname_author_year"
  # ## gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010", "asdasd") # must fail
  # gwas_ids <- c("CROHNS_Jostins2012", "LIPIDS_HDL_Teslovich2010") # must pass
  # return_as_df <- FALSE
  # 
  ### Check input
  allowed_styles <- c("fullname_author_year",
                      "fullname",
                      "abrv_author_year",
                      "abrv_year",
                      "abrv" # OBS: abrv may not be unique!
  )
  if (!style %in% allowed_styles) {
    stop("Got wrong style. Allowed styles are: [%s]", paste(allowed_styles, collapse=", "))
  }
  if (style == "abrv") {
    print("Warning: trait name abreviations are not unique.")
  }
  ### Read data
  df.gwas_database <- read_tsv(.data.rename_gwas)
  if (any(duplicated(df.gwas_database$gwas_id))) {
    stop("internal renaming data frame contains duplicates. Fix the .csv file")
  }
  ### Check that all gwas_ids are present in "data base"
  bool.found <- gwas_ids %in% df.gwas_database$gwas_id
  if (check_all_matches && !all(bool.found)) {
    gwas_ids_not_found <- gwas_ids[!bool.found]
    stop(sprintf("Some gwas_ids not found: [%s]", paste(gwas_ids_not_found, collapse=",")))
  }
  ### Format output conditional on style
  df.gwas_database.select <- df.gwas_database %>% filter(gwas_id %in% gwas_ids)
  if (style=="fullname_author_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_full, " (", author, ", ", year, ")"))
  } else if (style=="fullname") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_name_full))
  } else if (style=="abrv_author_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_abrv, " (", author, ", ", year, ")"))
  } else if (style=="abrv_year") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_abrv, " (", year, ")"))
  } else if (style=="abrv") {
    df.gwas_out_fmt <- df.gwas_database.select %>% 
      mutate(fmt=paste0(trait_abrv))
  } else {
    stop("Internal function error")
  }
  if (return_as_df) {
    gwas_out_fmt <- df.gwas_out_fmt
  } else(
    gwas_out_fmt <- df.gwas_out_fmt %>% pull(fmt) # returns vector
  )
  return(gwas_out_fmt)
}
