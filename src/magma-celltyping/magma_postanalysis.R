############### SYNOPSIS ###################
# Do MAGMA post analysis (e.g. Skene Fig2a 'GWAS heatmap')

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:



# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping"
setwd(wd)

# library(MAGMA.Celltyping) # loads EWCE package
library(tidyverse)


# ======================================================================= #
# ============================ CONSTANTS =============================== #
# ======================================================================= #

dataset <- "MACA"
# dataset <- "mousebrain"

# "out.magma_celltyping.XXX.results.MACA.csv"
# "out.magma_celltyping.SCZ_Ripke2014.results.mousebrain.csv"

dir.results <- "/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping"

# ======================================================================= #
# ============================ READ RESULTS =============================== #
# ======================================================================= #

list.files <- file.path(dir.results, list.files(path = dir.results, pattern = sprintf("*.results.%s.csv", dataset)))
list.dfs <- lapply(list.files, read_csv)
names(list.dfs) <- stringr::str_match(list.files, "out.magma_celltyping.(.*).results.*")[,2]
names(list.dfs)
df <- bind_rows(list.dfs, .id="gwas")
head(df)

# ======================================================================= #
# ============================ HEATMAP =============================== #
# ======================================================================= #


df.clean <- df %>% rename(annotation=COVAR, pval=P) %>%
  select(gwas, annotation, pval)
df.clean

df.clean <- df.clean %>% group_by(gwas) %>% mutate(
  rank = base::rank(pval), # rank() ranks in ascending order (smallest number == lowest rank [e.g. 1]). 
  pval_mlog10 = -log10(pval),
  significance = ifelse(pval < 0.05, "*", ""), # asterisk or empty string
  rank_txt = paste0(rank, significance)
) 
df.clean


p <- ggplot(df.clean, aes(x=gwas, y=annotation, fill=pval_mlog10)) + geom_tile() + geom_text(aes(label=rank_txt))
p <- p + labs(title=dataset) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + scale_fill_gradient2()
p
# ggsave(filename = "MAGMA_multigwas_heatmap.pdf", w=15, h=25)

