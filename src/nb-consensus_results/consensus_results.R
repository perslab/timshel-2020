############### SYNOPSIS ###################
# Combine results ('consus results') from multiple RolyPoly and MAGMA SEMs

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-consensus_results/"
setwd(wd)

# library(MAGMA.Celltyping) # loads EWCE package
library(tidyverse)


# ======================================================================= #
# ============================ CONSTANTS =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ READ MOUSEBRAIN METADATA =============================== #
# ======================================================================= #

file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/constants-annotation_name_mapping.mousebrain.csv"
df.metadata <- read_csv(file.metadata)
df.metadata

# ======================================================================= #
# ============================ MAGMA =============================== #
# ======================================================================= #

dir.results <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping/OUT"
gwas.name <- "BMI_Yengo2018"

list.files <- file.path(dir.results, list.files(path = dir.results, pattern = sprintf("(.*).%s.results.MOUSEBRAIN_EXT_(.*).csv", gwas.name)))
list.files
list.dfs <- lapply(list.files, read_csv)
names(list.dfs) <- stringr::str_match(list.files, "(.*).results.MOUSEBRAIN_EXT_(.*).csv")[,3] # ---> SEM model
names(list.dfs)
df <- bind_rows(list.dfs, .id="sem")
head(df)

### Spread (select columns first so that spread() will work)
df.spread <- df %>% select(sem, COVAR, P) %>% spread(key=sem, value=P)
df.spread

### Select specific SEMs, rank and calculate score
df.spread.rank <- df.spread %>% 
  select(COVAR, enrichment_log, z_score_pos_log, quantiles) %>% # select
  mutate_at(vars(-COVAR), funs(rank=rank)) %>%  # rank
  mutate(COVAR = df.spread$COVAR) %>% # re-adding COVAR
  select(COVAR, everything()) %>% #set COVAR first
  mutate(score=enrichment_log_rank+z_score_pos_log_rank+quantiles_rank) %>% # calculate score
  arrange(score)
df.spread.rank

### Copy
df.score.magma <- df.spread.rank

# ======================================================================= #
# ============================ ROLYPOLY  =============================== #
# ======================================================================= #

### RolyPoly - 'table.pvals.ALL_DATA.csv'
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/export-combined.rp_hm3.v3.nboot100/inference_rp_hm3.body_BMI_Yengo2018.tss.10kb.hm3.none.protein_coding_only.nboot100/table.pvals.ALL_DATA.csv"
df.rolypoly <- read_csv(file.in)
df.rolypoly %>% distinct(dataset)
df.rolypoly <- df.rolypoly %>% filter(dataset %in% c("mousebrain.enrichment_log", "mousebrain.z_score_pos_log", "mousebrain_skene.quantiles")) # filter 
df.rolypoly <- df.rolypoly %>% mutate(dataset=stringr::str_split_fixed(dataset, "\\.", 2)[,2])

### Rename some columns to be able to reuse the above MAGMA code [CONVENIENCE]
df <- df.rolypoly %>% rename(sem=dataset, COVAR=annotation, P=bp_value)

### COPY FROM ABOVE | Spread (select columns first so that spread() will work)
df.spread <- df %>% select(sem, COVAR, P) %>% spread(key=sem, value=P)
df.spread

### COPY FROM ABOVE | Select specific SEMs, rank and calculate score
df.spread.rank <- df.spread %>% 
  select(COVAR, enrichment_log, z_score_pos_log, quantiles) %>% # select
  mutate_at(vars(-COVAR), funs(rank=rank)) %>%  # rank
  mutate(COVAR = df.spread$COVAR) %>% # re-adding COVAR
  select(COVAR, everything()) %>% #set COVAR first
  mutate(score=enrichment_log_rank+z_score_pos_log_rank+quantiles_rank) %>% # calculate score
  arrange(score)
df.spread.rank

### Copy
df.score.rolypoly <- df.spread.rank

# ======================================================================= #
# ============================ COMBINE  =============================== #
# ======================================================================= #

### Join data frames
df.combined <- full_join(df.score.magma, df.score.rolypoly, by="COVAR", suffix=c(".magma", ".rolypoly"))
df.combined

### Add combined score
df.combined <- df.combined %>% mutate(score.combined = score.magma+score.rolypoly) %>% arrange(score.combined)

### Add metadata
df.combined <- full_join(df.combined, df.metadata, by=c("COVAR"="ClusterName"))
df.combined

### Write file
df.combined %>% write_csv("out.consensus_results.mousebrain.csv")

# ======================================================================= #
# ============================ PLOT  =============================== #
# ======================================================================= #

### Select and gather
df.combined.gather <- df.combined %>% select(COVAR, score.combined,
                              paste0(c("enrichment_log", "z_score_pos_log", "quantiles"), ".magma"), 
                              paste0(c("enrichment_log", "z_score_pos_log", "quantiles"), ".rolypoly")) %>% 
  gather(key="key", value="value", -COVAR, -score.combined)
df.combined.gather

### Order by score.combined
df.combined.gather <- df.combined.gather %>% arrange(score.combined) %>% mutate(COVAR = factor(COVAR, levels=unique(COVAR)))

### Plot
p <- ggplot(df.combined.gather, aes(x=COVAR, y=-log10(value), fill=score.combined)) + geom_col() + facet_wrap(~key, ncol=999) + coord_flip()
p
ggsave(filename="out.consensus_results.mousebrain.pdf", width=12, height=25)


