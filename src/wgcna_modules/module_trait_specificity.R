############### SYNOPSIS ###################
# Module trait specificity


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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

setwd(here("src/wgcna_modules"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

prefix_genomic_annot <- "wgcna.mousebrain-190213.fdr_sign_celltypes.continuous"

# ======================================================================= #
# ============================ FUNCTIONS =============================== #
# ======================================================================= #


load_ldsc_cts_results_wgcna <- function(file.ldsc_cts) {
  df.ldsc_cts <- read_tsv(file.ldsc_cts)
  mat_tmp_split_str <- stringr::str_split_fixed(df.ldsc_cts$Name,pattern="\\__",n=Inf)
  df.ldsc_cts <- df.ldsc_cts %>% mutate(module_id=mat_tmp_split_str[,ncol(mat_tmp_split_str)]) # add module ID from 'Name' string, e.g. maca_tissue_cell_type.slateblue4
  df.ldsc_cts <- df.ldsc_cts %>% rename(BETA=Coefficient, SE=Coefficient_std_error, P=Coefficient_P_value)
  return(df.ldsc_cts)
}

# ======================================================================= #
# ======================== LOAD LDSC CTS RESULTS ======================== #
# ======================================================================= #


### Load - MULTI GWAS
dir.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern=sprintf("%s.(.*).cell_type_results.txt", prefix_genomic_annot))
filenames <- filenames[!grepl(pattern="__CONDITIONAL__", filenames, perl=T)] # exclude any conditional results
list.dfs <- lapply(file.path(dir.data, filenames), load_ldsc_cts_results_wgcna)
names(list.dfs) <- stringr::str_match(filenames, pattern=sprintf("%s__(.*).cell_type_results.txt", prefix_genomic_annot))[,2] # ALT: filenames
names(list.dfs)
df.ldsc_cts <- list.dfs %>% bind_rows(.id="gwas")

# ======================================================================= #
# ================================ EXPORT =============================== #
# ======================================================================= #

df.ldsc_cts <- df.ldsc_cts %>% filter(module_id %in% c("lavenderblush", "lightpink3")) # filter modules

### format file [*not* compatible with multi-GWAS plotting]
df.ldsc_cts.export <- df.ldsc_cts %>% select(gwas, P, module_id) %>% spread(key=module_id, value=P) # format file
df.ldsc_cts.export <- df.ldsc_cts.export %>% arrange(lavenderblush) %>% mutate(gwas=factor(gwas, levels=gwas)) # order factor
file.out <- sprintf("out.trait_specificity.%s.csv", prefix_genomic_annot)
df.ldsc_cts.export %>% write_csv(file.out)

# ======================================================================= #
# ================================ PLOT =============================== #
# ======================================================================= #
order.gwas <- df.ldsc_cts %>% filter(module_id=="lavenderblush") %>% arrange(P) %>% pull(gwas)
df.plot <- df.ldsc_cts %>% select(gwas, P, module_id) %>% mutate(gwas = factor(gwas, levels=order.gwas))
df.plot <- df.plot %>% mutate(mlog10_P=-log10(P))
df.plot <- df.plot %>% filter(module_id=="lavenderblush")
df.plot <- df.plot %>% filter(!gwas %in% c("WHR_Pulit2019", "BMI_Pulit2019", "BMI_male_Pulit2019", "BMI_female_Pulit2019", "BMI_UPDATE_Yengo2018", "BMI_Locke2015"))

fdr_threshold <- 0.05/nrow(df.plot)
fdr_threshold

p <- ggplot(df.plot, aes(x=gwas, y=mlog10_P)) + 
  geom_col(position=position_dodge()) + 
  geom_col(data=df.plot %>% filter(P<fdr_threshold), fill="red", position=position_dodge()) + 
  geom_hline(yintercept=-log10(fdr_threshold), color="gray", linetype="dashed") +
  labs(x="", y=expression(-log[10](P))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.margin = unit(c(1,1,1,2), "cm"))
  # gghighlight(mlog10_P<-log10(0.05)) +
  # coord_flip()
p



