############### SYNOPSIS ###################
### AIM: XXX top ESmu genes for selected annotations

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

# library(knitr)
# library(kableExtra)
# library(docxtools)


source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ================================ LOAD DATA ============================ #
# ======================================================================= #

### Get ES data
load(here("out/es/mousebrain_all.es_obj.RData"))


### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

### Read and add meta data
df.metadata <- get_metadata(dataset_prefix="mousebrain_all")
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts

# ======================================================================= #
# ======================== EXTRACT AND PROCESS DATA ===================== #
# ======================================================================= #

### Filter: BMI results and prioritized cell-types
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
df <- df.ldsc_cts %>% 
  filter(gwas == "BMI_UKBB_Loh2018") %>%
  filter(annotation %in% filter.annotations)


cols_to_keep <- c("p.value",
                  "annotation",
                  "Description",
                  "NCells",
                  "Neurotransmitter",
                  "Probable_location",
                  "Region")
# POTENTIAL estimate,std.error,TaxonomyRank4
df <- df %>% select(cols_to_keep)

df <- df %>% mutate(Region = case_when(
  Region == "Midbrain dorsal" ~ "Midbrain",
  Region == "Midbrain dorsal,Midbrain ventral" ~ "Midbrain",
  Region == "Hippocampus,Cortex" ~ "Hippocampus/Cortex",
  TRUE ~ as.character(Region))
)

df <- df %>% mutate(Probable_location = str_replace_all(Probable_location, ", ", "\n"))
df <- df %>% mutate(Neurotransmitter = str_replace_all(Neurotransmitter, ", ", "\n"))

### docxtools: formatting numbers
### https://graphdr.github.io/docxtools/
### https://cran.r-project.org/web/packages/docxtools/vignettes/numbers-in-engineering-format.html
# df <- df %>% mutate("P-value"= df %>% select(p.value) %>% format_engr(.) %>% pull(p.value)) # # --> returns latex expression

# ======================================================================= #
# =========== MAIN TABLE: metadata + es genes for BMI cell-types ======== #
# ======================================================================= #

### Top ES genes
df.top_n <- get_annotation_es.top_n(sem_obj, annotations=filter.annotations, n_top_genes=10, es_metric="es_mu")
df.top_n <- df.top_n %>% select(annotation, es_genes_top_fmt) %>% distinct()
df.top_n

### Join with df
df <- df %>% left_join(df.top_n, by="annotation")

### Write out
df <- df %>% arrange(p.value)
file.out <- "tables/table-mousebrain_bmi_celltypes.csv"
df %>% write_csv(file.out)

# ======================================================================= #
# ===================== SOM table (csv): all genes ES =================== #
# ======================================================================= #

### ES table for all genes (SOM)
df.all_table <- get_annotation_es.table(sem_obj, annotations=filter.annotations, es_metric="es_mu")
file.out <- "tables/table-es_mu.mousebrain_bmi_celltypes.csv.gz"
df.all_table %>% write_csv(file.out)


# ======================================================================= #
# ================================ LEFTOVER FUNCTIONS ============================ #
# ======================================================================= #

scientific_10x <- function(values, digits = 1) {
  # COPY FROM https://rdrr.io/github/jpshanno/ecoflux/src/R/Helper_Functions.R
  # ---> return expression()
  if(!is.numeric(values)){
    stop("values must be numbers")
  }
  if(grepl("^\\d{2}$", digits)){
    stop("digits must a one or two digit whole number")
  }
  
  x <- sprintf(paste0("%.", digits, "e"), values)
  
  x <- gsub("^(.*)e", "'\\1'e", x)
  
  longestExponent <- max(sapply(gregexpr("\\d{1,}$", x), attr, 'match.length'))
  zeroTrimmed <- ifelse(longestExponent > 2,
                        paste0("\\1", paste(rep("~", times = longestExponent-1), collapse = "")),
                        "\\1")
  x <- gsub("(e[+|-])[0]", zeroTrimmed, x)
  
  x <- gsub("e", "~x~10^", x)
  
  if(any(grepl("\\^\\-", x))){
    x <- gsub("\\^\\+", "\\^~~", x)
  } else {
    x <- gsub("\\^\\+", "\\^", x)
  }
  # return this as an expression
  parse(text=x)
} 

# scientific_10x(df$p.value)
# df %>% mutate("P-value"=sapply(df$p.value,scientific_10x,simplify=F))
# df %>% mutate("P-value"=unlist(lapply(df$p.value,scientific_10x)))


