############### SYNOPSIS ###################
### AIM: plot mousebrain meta-data


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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


setwd(here("src/GE-mousebrain"))


# ======================================================================= #
# =========================== READ METADATA =========================== #
# ======================================================================= #

### Read
file.metadata <- here("data/expression/mousebrain/mousebrain-agg_L5.metadata.csv")
df.metadata <- read_csv(file.metadata)

# ======================================================================= #
# =========================== PLOT NT METADATA ========================== #
# ======================================================================= #

### Wrangle
df.metadata.nt <- df.metadata %>% 
  select(annotation, starts_with("nt.")) %>%
  gather(key="nt", value="value", -annotation) %>%
  mutate(nt=str_replace_all(nt, "nt.", ""))


### PLOT
ggplot(df.metadata.nt, aes(x=annotation, y=nt, color=nt)) + 
  geom_point(data=df.metadata.nt %>% filter(value==1)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle=90))

