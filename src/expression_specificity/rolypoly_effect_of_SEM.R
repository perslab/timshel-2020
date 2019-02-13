
library(tidyverse)



# =========================================================================================================== #
# ======================= EXPLPORE PROPERTIES OF ROLYPOLY *different 'SEM'* MOUSEBRAIN  RESULTS ================ #
# =========================================================================================================== #

rm(list=ls()) # clear workspace

library(GGally) # for ggpairs and ggcorr
library(ggridges)

### Load 'table.pvals.ALL_DATA.csv'
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/RP-meta/export-combined.rp_hm3.v3.nboot100/inference_rp_hm3.body_BMI_Yengo2018.tss.10kb.hm3.none.protein_coding_only.nboot100/table.pvals.ALL_DATA.csv"
df <- read_csv(file.in)
df <- df %>% filter(grepl("mousebrain",dataset))
df <- df %>% mutate(dataset=stringr::str_replace(dataset, ".*\\.", "")) # shorten dataset names (rename)
# ALTERNATIVE SHORTEN? del: lapply(stringr::str_split(colnames(df.clean), "\\."), '[[', -1)
# ALTERNATIVE SHORTEN? del: stringr::str_match(colnames(df.clean), "\\.?(.*?)$") # match

# ======================= Matrix plot of beta ======================= #
df.beta <- df %>% select(annotation, dataset, bootstrap_estimate)  %>% spread(key=dataset, value=bootstrap_estimate)
### Matrix plot
p <- df.beta %>% select(-annotation) %>% ggpairs()
p


# ======================= Ridge plot of SE / beta ======================= #
df.se <- df %>% select(annotation, dataset, bootstrap_error)
df.se
p <- ggplot(df.se, aes(x=bootstrap_error, y=dataset)) + geom_density_ridges()
p
p + xlim(-0.0000005, 1e-5)

### beta + se combined
df.beta_se <- df %>% select(annotation, dataset, bootstrap_error, bootstrap_estimate) %>% gather(key="statistic", value="bootstrap_error_and_beta", bootstrap_estimate, bootstrap_error)
df.beta_se
p <- ggplot(df.beta_se, aes(x=bootstrap_error_and_beta, y=dataset, fill=statistic)) + geom_density_ridges() #+ facet_wrap(~statistic)
p + xlim(-1e-6, 1e-5)

df.t_value <- df %>% select(annotation, dataset, bt_value)
df.t_value
p <- ggplot(df.t_value, aes(x=bt_value, y=dataset)) + geom_vline(xintercept=0) + geom_density_ridges()
p
p + xlim(-1e-6, 1e-5)



# =========================================================================================================== #
# =======================  TMP TMP TMP TMP TMP ================ #
# =========================================================================================================== #

# ======================= TMP CLEAN INFERENCE - MAY DELETE ======================= #
files <- c("body_BMI_Yengo2018",
  "lipids_LDL_Willer2013",
  "mental_SCZ_Ripke2014")
for (f in files) {
  # f <- "body_BMI_Yengo2018" # test
  # file.in <- sprintf("/scratch/sc-genetics/out.rolypoly_objs-v3.hm3_protein_coding_only.univariate-nboot100/rolypoly_objs.%s.tss.10kb.hm3.none.protein_coding_only.nboot100.inference.RData", f)
  print(file.in)
  # load(file.in)
  names(list.rp_inference)
  list.rp_inference[["mousebrain.enrichment_quantiles"]] <- NULL
  list.rp_inference[["mousebrain_skene.specificity_quantiles"]] <- NULL
  list.rp_inference[["mousebrain.z_score_quantiles"]] <- NULL
  # base::save(list.rp_inference, list.run_parameters,file=file.in)
}
