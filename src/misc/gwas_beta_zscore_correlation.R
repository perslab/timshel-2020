

### cov_SMOKING_STATUS (example of a binary trait) : checking if all betas are > 0, meaning OR
df.alkes_ukbb <- read_tsv("/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/cov_SMOKING_STATUS.sumstats.gz")
names(df.alkes_ukbb)
summary(df.alkes_ukbb$Beta) # --> 
# M in.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.4391500 -0.0030290 -0.0000188 -0.0000337  0.0029946  0.3123350 

### LDSC format | LOG TRANS BETA with some ZERO vals
df <- read_tsv("/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/alkesgroup/sumstats_formatted/Ever_Smoked.sumstats")
summary(df$Z) # min = 0
anyNA(df$Z)
anyNA(log(df$Z)) # log(0) ==> Inf
df %>% filter(df$Z==0)


#### ROLYPOLY FMT | R models NA values as character
df.fmt <- read_tsv("/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/alkesgroup-snpsnap_mapped/Ever_Smoked.gwassumstats.rolypoly_fmt.tab.gz")
class(df.fmt$se[1])
x<-df.fmt$se[1]
is.na(df.fmt$se[1])



# ======================= GWAS beta vs Z correlation ======================= #

df.gwas <- read_tsv("/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/SWB_Full.txt.gz")
plot(df.gwas$Beta, -log10(df.gwas$Pval))
     

### Edu
df.edu.raw <- read_tsv("/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA2_Okbay2016/EduYears_Main.txt.gz")
df.edu.munge <- read_tsv("/raid5/projects/timshel/sc-genetics/ldsc/src/tutorial-james/EA2_excl_23andMe_noGC.sumstats.gz")

df.edu <- left_join(df.edu.munge, df.edu.raw, by=c("SNP"="MarkerName"))

cor(df.edu$Z, df.edu$Beta, use="complete.obs") # ---> 0.9052704
# plot(df.edu$Z, df.edu$Beta)

cor(abs(df.edu$Z), -log10(df.edu$Pval), use="complete.obs") # ---> 0.9442885
# plot(abs(df.edu$Z), -log10(df.edu$Pval))



