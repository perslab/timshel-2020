############### SYNOPSIS ###################
# Plot WGCNA module metadata correlation for Novo bulk analysis.

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)

# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-novo_bulk/"
setwd(wd)

# ======================================================================= #
# ============================ PLOT 1: RHO BAR PLOT ============================== #
# ======================================================================= #

dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables"

### read expression data - NOVO BULK | rho (you can also do the same for log_fdr)
filenames <- list.files(path=dir.data,  pattern="nn_lira_sema_per_brain_area_run1_.*_metadata_corr_rho.csv") # RHO
# ^ match | nn_lira_sema_per_brain_area_run1_NTS_metadata_corr_rho.csv
filenames
list.dfs <- lapply(file.path(dir.data, filenames), read_csv)
names(list.dfs) <- filenames # e.g. xxx
names(list.dfs)
df <- list.dfs %>% purrr::reduce(full_join, by = "X1")
# df.pvals <- read_csv("/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_all_metadata_corr_logfdr.csv")
# df.pvals
df <- df %>% rename(treatment = X1)
df <- df %>% filter(treatment %in% c("Vehicle",
                        "Liraglutide",
                        "Semaglutide",
                        "WM"))
df.spread <- df %>% gather(key="key", value="rho", -treatment) 
df.spread

### plot
p <- ggplot(df.spread, aes(x=key, y=rho, fill=treatment)) + geom_col(position = position_dodge(0.3))
p <- p + theme(axis.text.x=element_text(angle = 45, hjust = 1))
p

# ggsave(filename = "novo_bulk_wgcna_correlation_plot.pdf", width = 15, height = 8)


# ======================================================================= #
# ============================ PLOT 2 - with P-val significance ============================== #
# ======================================================================= #

dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables"

### read rho (you can also do the same for log_fdr)
filenames <- list.files(path=dir.data,  pattern="nn_lira_sema_per_brain_area_run1_.*_metadata_corr_rho.csv") # RHO
# ^ match | nn_lira_sema_per_brain_area_run1_NTS_metadata_corr_rho.csv
filenames
list.dfs <- lapply(file.path(dir.data, filenames), read_csv)
names(list.dfs) <- filenames # e.g. xxx
names(list.dfs)
df.rho <- list.dfs %>% purrr::reduce(full_join, by = "X1")
### read P-vals
df.pval <- read_csv("/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_all_metadata_corr_logfdr.csv")
df.pval
### join
df <- full_join(df.rho, df.pval, by="X1", suffix = c("-DUMMY_JOIN_STRING-rho", "-DUMMY_JOIN_STRING-pval"))
### filter
df <- df %>% rename(treatment = X1)
df <- df %>% filter(treatment %in% c("Vehicle",
                                     "Liraglutide",
                                     "Semaglutide",
                                     "WM"))
df.gather1 <- df %>% gather(key="key", value="value", -treatment) 
df.gather2 <- df.gather1 %>% separate(col="key", into=c("module", "key_stat"), sep="-DUMMY_JOIN_STRING-") # split 'key' column by DUMMY JOIN STRING
df.gather3 <- df.gather2 %>% spread(key="key_stat", value="value")
df.gather3

### plot x-y correlation between pval and rho
ggplot(df.gather3, aes(x=pval, y=rho)) + geom_point()

### plot
p <- ggplot(df.gather3, aes(x=module, y=rho, fill=treatment)) + 
  geom_col(position = position_dodge(0.3)) +
  geom_point(aes(color=treatment, size=pval)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
p



# ======================================================================= #
# ============================ XXXXXXXXXX ============================== #
# ======================================================================= #



# ======================================================================= #
# ============================ LEFTOVERS ============================== #
# ======================================================================= #


### Alternative READING DATA | using bind_cols, but I hard problmes making it give the desired output.
# list.dfs <- lapply(file.path(dir.data, filenames), read_csv)
# list.dfs <- lapply(list.dfs, function(df) {df %>% select(-X1)}) # remove first column, otherwise bind_cols() result will be 'messy'
# names(list.dfs) <- filenames # e.g. xxx
# names(list.dfs)
# df <- bind_cols(list.dfs) # BINDS PER POSITION, so all data.frames must have rows in the same order
# df

# filenames_shorten <- stringr::str_match(filenames, "nn_lira_sema_per_brain_area_run1_(.*)_metadata_corr_(.*).csv")[,3] # e.g. 
# filenames_shorten
# ^ match | nn_lira_sema_per_brain_area_run1_NTS_metadata_corr_logfdr.csv
