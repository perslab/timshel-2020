############### SYNOPSIS ###################
### AIM: Parse mousebrain metadata.
# 1: neurotransmitters
# 2: clean Region

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

setwd(here("src/datasets-expression/mousebrain"))

# ======================================================================= #
# ================================ READ DATA ================================ #
# ======================================================================= #

cols_metadata_keep <- c("ClusterName",
                        "Class",
                        "Description",
                        "NCells",
                        "Neurotransmitter",
                        "Probable_location",
                        "Region",
                        "LeafOrder", # order in taxonomy Mousebrain Fig1C (the taxonomy plot). Same order as in Table S3
                        "TaxonomyRank1",
                        "TaxonomyRank2",
                        "TaxonomyRank3",
                        "TaxonomyRank4",
                        "Comment"
                        )

### Read
file.metadata <- here("data/expression/mousebrain/mousebrain_agg_L5.metadata_full.csv") # exported via loompy from "L5_All.agg.loom"
df.metadata <- read_csv(file.metadata)

### Select metadata
df.metadata <- df.metadata %>% 
  select(cols_metadata_keep) %>% 
  rename(annotation = ClusterName)

# ======================================================================= #
# ======================= ADD TAXONOMY METADATA ========================= #
# ======================================================================= #
# We have created additional taxonomy metadata for mousebrain data
# It consists of a mapping from "TaxonomyRank4" to "TaxonomyRank4_reduced1" and "tax_order_idx_mb_fig1c".
  # TaxonomyRank4_reduced1 (n_tax=24) is a more compact TaxonomyRank4 (n_tax)=39) more suitable for our analysis.
  # tax_order_idx_mb_fig1c orders the TaxonomyRank4_reduced1 by the Mousebrain publication Fig1C (the taxonomy plot)

### Read annotation taxonomy metadata
file.taxonomy_metadata <- here("data/expression/mousebrain/mousebrain_table_s5-taxonomy.csv")
df.taxonomy_metadata <- read_csv(file.taxonomy_metadata)


df.taxonomy_metadata <- df.taxonomy_metadata %>% select(TaxonomyRank4, TaxonomyRank4_reduced1, tax_order_idx_mb_fig1c)
if ( !(n_distinct(df.taxonomy_metadata$TaxonomyRank4_reduced1) == n_distinct(df.taxonomy_metadata$tax_order_idx_mb_fig1c)) ) {
  stop("Columns 'tax_order_idx_mb_fig1c' and 'TaxonomyRank4_reduced1' must have the same number of unique values. Otherwise downstream sorting of factors will fail.")
}


### Add meta data
df.metadata <- df.metadata %>% left_join(df.taxonomy_metadata, by="TaxonomyRank4") # add taxonomy meta data [*OBS*: joining by TaxonomyRank4 and not annotation]

# ======================================================================= #
# ============================ SELECTION ============================= #
# ======================================================================= #

df.metadata.select <- df.metadata 

# We here subset on neurons only because it makes it easier to interpret the counts from the Neurotransmitter classes when we compare to Linnarson assignments
### Subset on Neurons only
# df.metadata.select <- df.metadata %>% filter(Class=="Neurons")

### Replace NA values in Neurotransmitter. This makes everything downstream easier.
# df.metadata.select <- df.metadata.select %>% mutate(Neurotransmitter = if_else(is.na(Neurotransmitter), "Missing", Neurotransmitter))

# ======================================================================= #
# ============================ ASSIGN EX/IN ============================= #
# ======================================================================= #
# Here we assign neurotransmitters to Ex/In
# metabotropic and ionotropic receptors

class.ex <- c("Glutamate")
class.in <- c("GABA", "Glycine")
class.monoamines <- c("Adrenaline", 
                      "Noradrenaline", 
                      "Dopamine",
                      "Serotonin")
class.ach <- c("Acetylcholine")
class.no <- c("Nitric oxide")
# class.other <- c("Nitric oxide", "Acetylcholine")

### Make assignment
# *IMPORTANT*: because of how case_when works, the below order list the 'precedence' of transmitter class assignment
# "Like an if statement, the arguments to case_when are evaluated in order, so you must proceed from the most specific to the most general."
df.metadata.select <- df.metadata.select %>% 
  mutate(Neurotransmitter_class=case_when(
    grepl(pattern=paste(class.ex, collapse="|"), Neurotransmitter) & grepl(pattern=paste(class.in, collapse="|"), Neurotransmitter) ~  "Ex/In",
    grepl(pattern=paste(class.ex, collapse="|"), Neurotransmitter) ~ "Excitatory",
    grepl(pattern=paste(class.in, collapse="|"), Neurotransmitter) ~ "Inhibitory",
    grepl(pattern=paste(class.monoamines, collapse="|"), Neurotransmitter) ~ "Monoamines",
    grepl(pattern=paste(class.ach, collapse="|"), Neurotransmitter) ~ "Acetylcholine",
    grepl(pattern=paste(class.no, collapse="|"), Neurotransmitter) ~ "Nitric oxide",
    Class!="Neurons" ~ "Non-neuron", # this should come BEFORE is.na(Neurotransmitter)
    is.na(Neurotransmitter) ~ "Missing neurotransmitter", # 
    TRUE ~ "Other neurotransmitter" # In case our class.* vectors does not capture the transmitter
  )
  )

### Add Linnarsons (incomplete) asssignment based on "Description"
df.metadata.select <- df.metadata.select %>% 
  mutate(Neurotransmitter_class_linnarson=case_when(
    str_detect(Description, pattern=regex("Excitatory", ignore_case=T)) ~ "Ex",
    str_detect(Description, pattern=regex("Inhibitory", ignore_case=T)) ~ "In",
    # ^ str_detect() returns a logical vector.
    TRUE ~ NA_character_
  )
  )

### Cross validate
df.metadata.select %>% count(Neurotransmitter_class)
# 1 Acetylcholine                7
# 2 Ex/In                        4
# 3 Excitatory                 102
# 4 Inhibitory                  78
# 5 Missing neurotransmitter     5
# 6 Monoamines                  12
# 7 Nitric oxide                 6
# 8 Non-neuron                  51
df.metadata.select %>% count(Neurotransmitter_class_linnarson)
# 1 NA                                 106
# 2 Ex                                  60
# 3 In                                  48
df.metadata.select %>% count(Neurotransmitter_class, Neurotransmitter_class_linnarson)
# Neurotransmitter_class Neurotransmitter_class_linnarson     n
# 1 Ex                     NA                                  42
# 2 Ex                     Ex                                  60 OK
# 3 Ex/In                  NA                                   4
# 4 In                     NA                                  30
# 5 In                     In                                  48 OK
# 6 Monoamines             NA                                  12
# 7 Other                  NA                                  13
# 8 Undefined class        NA                                   5


# ======================================================================= #
# ==================== PRE-PROCESS Neurotransmitter ===================== #
# ======================================================================= #

### STEP 1 | Separate rows
df.nt <- df.metadata.select %>% 
  mutate(Neurotransmitter_split = str_split(Neurotransmitter, pattern=regex(",\\s*(?![^()]*\\))"))) %>% # split by comma, but not inside parentheses. REF: https://stackoverflow.com/a/26634150
  # ^ we cannot use separate_rows() because it does not accept regex. # REF: https://stackoverflow.com/a/31514711
  unnest(Neurotransmitter_split) # unnest our list column
# mutate(Neurotransmitter = str_trim(Neurotransmitter, side="both")) # trim any leading/trailing whitespace.
df.nt

### OPTIONAL: Count
df.nt.count <- df.nt %>% count(Neurotransmitter_split) %>% arrange(desc(n))
df.nt.count
# Neurotransmitter_split         n
# <chr>                      <int>
# 1 GABA                          82
# 2 Nitric oxide                  47
# 3 Glutamate (VGLUT2)            42
# 4 Glutamate (VGLUT1)            36
# 5 Acetylcholine                 17
# 6 Glutamate (VGLUT1,VGLUT2)     13
# 7 Glutamate (VGLUT3)            11
# 8 Noradrenaline                  8
# 9 Glycine (GLYT1,GLYT2)          6
# 10 Glycine (GLYT2)                6
# 11 NA                             5
# 12 Serotonin                      5
# 13 Dopamine                       4
# 14 Glycine (GLYT1)                2
# 15 Adrenaline                     1
# 16 Glutamate (VGLUT1, VGLUT3)     1
# 17 Glutamate (VGLUT1,VGLUT3)      1
# 18 Glutamate (VGLUT2, VGLUT3)     1
# 19 Glutamate (VGLUT2,VGLUT3)      1

### STEP 2 | Extract transmitter information inside parentheses (needed for e.g. "Glutamate (VGLUT1, VGLUT3)" or "Glutamate (VGLUT1)")
df.nt <- df.nt %>% mutate(Neurotransmitter_split_extract = str_match(Neurotransmitter_split, pattern="\\((.*)\\)")[,2])
# ^ Neurotransmitter_split_extract will have NA values for entries that does NOT contain parentheses (e.g. "Adrenaline" or "Nitric oxide)
df.nt <- df.nt %>% mutate(Neurotransmitter_split_extract = if_else(is.na(Neurotransmitter_split_extract), Neurotransmitter_split, Neurotransmitter_split_extract)) # If NA then assign Neurotransmitter_split value
# ^ Here we 'merge' Neurotransmitter_split_extract and Neurotransmitter_split.
### Separate rows: e.g. split "VGLUT1,VGLUT2" into two rows.
df.nt <- df.nt %>% 
  separate_rows(Neurotransmitter_split_extract, sep = ",") %>% 
  mutate(Neurotransmitter_split_extract = str_trim(Neurotransmitter_split_extract, side="both")) # trim any leading/trailing whitespace

head(df.nt %>% select(annotation, Neurotransmitter_split_extract))

### STEP3 | Spread: to make a [annotation x neurotransmitter] data frame.
df.nt.spread <- df.nt %>% 
  select(annotation, Neurotransmitter_split_extract) %>%
  mutate(value=1) %>% # this is needed as a 'dummy' value
  spread(key=Neurotransmitter_split_extract, value=value)
df.nt.spread
# annotation Acetylcholine Adrenaline Dopamine  GABA GLYT1 GLYT2 Missing `Nitric oxide` Noradrenaline Serotonin VGLUT1 VGLUT2 VGLUT3
# 1 CBGRC                 NA         NA       NA    NA    NA    NA      NA             NA            NA        NA      1     NA     NA
# 2 CBINH1                NA         NA       NA     1    NA    NA      NA              1            NA        NA     NA     NA     NA
### Add Neurotransmitter_class
# df.nt.spread <- df.nt.spread %>% left_join(df.metadata.select %>% select(annotation, Neurotransmitter_class, Neurotransmitter_class_linnarson), by="annotation")
### Add prefix to neurotransmitter columns.
df.nt.spread <- df.nt.spread %>% rename_at(vars(-annotation), ~ paste0("nt.", .))
df.nt.spread


### STEP 4 | Add Neurotransmitter_class info by joining df.nt.spread ('matrix') and df.metadata.select (EX/IN assignment).
df.metadata.select <- df.metadata.select %>% left_join(df.nt.spread, by="annotation")
# ^ this is currently not used because we export all annotations (df.metadata)

# ======================================================================= #
# ============================== REGION ================================= #
# ======================================================================= #

### Clean Names
df.metadata <- df.metadata %>% 
  mutate(Region = str_replace_all(Region, pattern="Medullae", replacement="Medulla"))
# ^ One cell-type (CBNBL1) was annotated with "Medullae"

# ======================================================================= #
# ============================== EXPORT ================================= #
# ======================================================================= #

### Join with all annotations
df.metadata <- df.metadata %>% left_join(df.metadata.select %>% select(annotation, Neurotransmitter_class, Neurotransmitter_class_linnarson), by="annotation") # add Neurotransmitter_class columns
df.metadata <- df.metadata %>% left_join(df.nt.spread, by="annotation") # add nt.* columns

### Write file
file.metadata_out <- here("data/expression/mousebrain/mousebrain.metadata.csv")
df.metadata %>% write_csv(file.metadata_out)
