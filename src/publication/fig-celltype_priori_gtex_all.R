############### SYNOPSIS ###################
### AIM: Make GTEx prioritzation bar plot ('DEPICT style')

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

library(patchwork)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.results)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id %in% "gtex_v8_filtered_normLog_all")
df.ldsc_cts %>% count(specificity_id)

# ======================================================================= #
# ============================ BMI point plot =============================== #
# ======================================================================= #

### Create 'plotting' data frame
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results

### Add fdr_significant flag (within GWAS)
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))

### Add columns
# df.ldsc_cts <- df.ldsc_cts %>% mutate(p.value.mlog10 = -log10(p.value))

# ======================================================================= #
# ============================ Add meta data =============================== #
# ======================================================================= #

df.ldsc_cts <- df.ldsc_cts %>% left_join(get_metadata("gtex_all"), by="annotation")
df.ldsc_cts


# ======================================================================= #
# ====================== FORMAT TO GTEx 'DEPICT fig' ====================== #
# ======================================================================= #

# "Name","Nominal.P.value","False.discovery.rate

df.input.raw <- df.ldsc_cts %>% 
  mutate(tmp=if_else(fdr_significant, "<0.01", ">=0.20"),
         Label=str_replace_all(tissue_sub, "^(.*) - ", "")) %>%
  select(
  "Name"=tissue_sub,
  "Nominal.P.value"=p.value,
  "False.discovery.rate"=tmp,
  Label
)

df.input.raw


# ======================================================================= #
# =============================== FUNCTION CALL ========================= #
# ======================================================================= #

source(here("src/publication/libx-gtex_figure.R"))


df.tissue_enrichment <- function.input_procces(df=as.data.frame(df.input.raw), cols2keep=colnames.gtex)
df.tissue_enrichment
### Updating df
df.tissue_enrichment <- function.gtex.set_variables(df.tissue_enrichment)
str(df.tissue_enrichment)

### Plot
p.gtex <- function.plot.barplot(df.tissue_enrichment, xlabel="Tissues")
p.gtex

### Adjust
p.gtex <- p.gtex + theme(plot.margin = unit(c(2,1,1,1), "cm")) # (t, r, b, l) widen margin
p.gtex

### export
file.out <- "figs/fig_celltypepriori.gtex.pdf"
ggsave(plot=p.gtex, filename=file.out, width=8, height=4)



# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #
