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

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/wgcna_modules"))



# ======================================================================= #
# ================================ READ DATA =============================== #
# ======================================================================= #

file.data <- here("results/prioritization_modules--fdr_sign_celltypes.multi_gwas.csv.gz")
df.ldsc_cts <-  read_csv(file.data)


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



