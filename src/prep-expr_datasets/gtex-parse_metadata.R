############### SYNOPSIS ###################
# Parse GTEx meta data


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ============================== SETUP =============================== #
# ======================================================================= #

library(tidyverse)
library(here)

# ======================================================================= #
# ============================== CONSTANTS =============================== #
# ======================================================================= #


if (!dir.exists(here("tmp-data"))) dir.create(here("tmp-data"))
if (!dir.exists(here("tmp-data","expression"))) dir.create(here("tmp-data", "expression"))
dir.data <- here("tmp-data/expression")

# ======================================================================= #
# ============================== DOWNLOAD =============================== #
# ======================================================================= #

downloadURLraw <- "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
download.file(downloadURLraw, destfile=file.path(dir.data, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))

# ======================================================================= #
# ======================== PARSE METDA DATA ======================= #
# ======================================================================= #

### READ Sample Annotations
file.annotations <- file.path(dir.data, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
df.annot <- read_delim(file.annotations, delim="\t")
problems(df.annot) # don't worry about problems with SMGTC column. It does not matter

df.annot <- df.annot %>% select(sample_id=SAMPID, tissue=SMTS, tissue_sub=SMTSD, RIN=SMRIN)
df.annot %>% count(tissue) %>% arrange(desc(n))
df.annot %>% count(tissue_sub) %>% arrange(desc(n))

### Make meta data
df.metadata <- df.annot %>% distinct(tissue, tissue_sub)
df.metadata <- df.metadata %>% mutate(annotation=str_replace_all(tolower(make.names(tissue_sub)), "(\\.)+", "_"))
df.metadata <- df.metadata %>% mutate(annotation=str_replace_all(annotation, "_$", ""))
df.metadata <- df.metadata %>% select(annotation, everything())

### Export
df.metadata %>% write_csv(here("data/expression/gtex/gtex-metadata.csv"))
