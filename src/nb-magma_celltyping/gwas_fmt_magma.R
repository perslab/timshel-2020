############### SYNOPSIS ###################
# REFORMAT GWAS TO MAGMA

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping"
setwd(wd)

### Output file
# SNP     P       N

DIR_OUT <- "/scratch/tmp-magma_gwas"




# ======================================================================= #
# ============================  UKBB GWAS (e.g. BLOOD) =============================== #
# ======================================================================= #


### blood_EOSINOPHIL_COUNT
file.out <- file.path(DIR_OUT, "blood_EOSINOPHIL_COUNT.txt")
N_study <- 444656
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/blood_EOSINOPHIL_COUNT.sumstats.gz"
df <- read_tsv(file.in)
#df <- df %>% rename(SNP=rsID, P=pval)
#df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df <- df %>% filter(complete.cases(.)) # OBS: remove NA values
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

### SNIPPET
# SNP     CHR     POS     A1      A2      REF     EAF     Beta    se      P       N       INFO
# rs10399793      1       49298   T       C       T       0.376299        -0.00793522     0.00336227      2.2E-02 444656  0.342797
# rs2462492       1       54676   C       T       C       0.599381        -0.0025273      0.00333118      5.3E-01 444656  0.340158
# rs3107975       1       55326   T       C       T       0.991547        0.0113799       0.0185189       4.1E-01 444656  0.324228



# ======================================================================= #
# ============================  ROSA GWAS =============================== #
# ======================================================================= #

### Daytime_dozing_sleeping (contains X-chromosome variants, but they will be removed)
file.out <- file.path(DIR_OUT, "Daytime_dozing_sleeping_UKBB2018.txt")
N_study <- 360855
file.in <- "/projects/rosa/gwas_corr/gwas_data/Daytime_dozing_sleeping_1220.tab.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=rsID, P=pval)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df <- df %>% filter(complete.cases(.)) # OBS: remove NA values
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

# Nap_during_day (contains X-chromosome variants, but they will be removed)
file.out <- file.path(DIR_OUT, "Nap_during_day_UKBB2018.txt")
N_study <- 359752
file.in <- "/projects/rosa/gwas_corr/gwas_data/Nap_during_day_1190.tab.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=rsID, P=pval)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df <- df %>% filter(complete.cases(.)) # OBS: remove NA values
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

# rsID    beta    se      pval    snp_maf chr     pos
# rs547522712     197.35  124.557 0.113101        5.45038e-09     1       15791
# rs568226429     -0.291333       0.243978        0.232443        5.77195e-06     1       69487
# rs2531267       -0.0783097      0.0437563       0.0735062       0.000188474     1       69569

### > problems(df) ---> DATA CONTAIN X-chromosome variants. This is not a problem: SNPs will get NA in 'chr' column and will be removed
# # A tibble: 427,164 x 5
# row col   expected   actual file                                                                    
# <int> <chr> <chr>      <chr>  <chr>                                                                   
#   1 13364304 chr   an integer X      '/projects/rosa/gwas_corr/gwas_data/Daytime_dozing_sleeping_1220.tab.gz'
# 2 13364305 chr   an integer X      '/projects/rosa/gwas_corr/gwas_data/Daytime_dozing_sleeping_1220.tab.gz'
# 3 13364306 chr   an integer X      '/projects/rosa/gwas_corr/gwas_data/Daytime_dozing_sleeping_1220.tab.gz'
# 4 13364307 chr   an integer X      '/projects/rosa/gwas_corr/gwas_data/Daytime_dozing_sleeping_1220.tab.gz'
# 5 13364308 chr   an integer X      '/projects/rosa/gwas_corr/gwas_data/Daytime_dozing_sleeping_1220.tab.gz'

# ======================================================================= #
# ============================  NULL GWAS  =============================== #
# ======================================================================= #
### TEST
# head /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NULL_GWAS/plink_out.1KG_phase3_EUR_null_gwas.P1.qassoc | perl -lane 'if ($.==1) {print "SNP\tP\tN"; next}; print "$F[1]\t$F[8]\t503"'

# N_study = 503 (individuals from 1KGP)

### BASH command ---> WORKS
# for i in {1..10}
# do
# cat /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NULL_GWAS/plink_out.1KG_phase3_EUR_null_gwas.P${i}.qassoc | perl -lane 'if ($.==1) {print "SNP\tP\tN"; next}; print "$F[1]\t$F[8]\t503"' > /scratch/tmp-magma_gwas/1KG_phase3_EUR_null_gwas_P${i}.txt
# done


### DOES NOT WORK (needs more tweeking with the quotes but not worth it)
# null_gwas_idx <- 1:10
# for (i in null_gwas_idx) {
#   perl_cmd <- 'if ($.==1) {print "SNP\tP\tN"; next}; print "$F[1]\t$F[8]\t503"'
#   cmd <- sprintf("head /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NULL_GWAS/plink_out.1KG_phase3_EUR_null_gwas.P%s.qassoc | perl -lane '%s' > /scratch/tmp-magma_gwas/1KG_phase3_EUR_null_gwas.P%s.txt", i, perl_cmd, i)
#   system2(cmd)
#   print(cmd)
# }


# ======================================================================= #
# ============================  FI_FG_Lagou2018  =============================== #
# ======================================================================= #


# snp     effect_allele   other_allele    eaf_hapmap_CEU  male_beta       male_se male_pvalue     female_beta     female_se       female_pvalue   sex_diff_pvalue sex_het_pvalue
# rs10000010      C       T       .446    -3.1e-03        3.0e-03 .302412 4.6e-04 2.8e-03 .872612 .580085 .389957
# rs10000012      C       G       .854    7.9e-04 5.7e-03 .889815 1.4e-02 5.1e-03 .005137 .019697 .078337
# rs1000003       A       G       .872    5.1e-03 5.2e-03 .33453  2.8e-03 4.7e-03 .555295 .527504 .750353
# rs10000037      A       G       .257    1.0e+00 4.5e-03 .915173 1.0e+00 4.1e-03 .978019 .993969 .951392


### FG | MALE
# fasting glucose (up to 67,506 men/73,089 women) 
file.out <- file.path(DIR_OUT, "FG_male_Lagou2018.txt")
N_study <- 67506 # MALE SAMPLE SIZE
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FG_STAGE1_2_3_SEX_GWAS_2018.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=snp, P=male_pvalue)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

### FI | MALE
# fasting insulin (up to 47,806 men/50,404 women),
file.out <- file.path(DIR_OUT, "FI_male_Lagou2018.txt")
N_study <- 47806 # MALE SAMPLE SIZE
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FI_STAGE1_2_3_SEX_GWAS_2018.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=snp, P=male_pvalue)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)



### FG | FEMALE
# fasting glucose (up to 67,506 men/73,089 women) 
file.out <- file.path(DIR_OUT, "FG_female_Lagou2018.txt")
N_study <- 73089 # FEMALE SAMPLE SIZE
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FG_STAGE1_2_3_SEX_GWAS_2018.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=snp, P=female_pvalue)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

### FI | FEMALE
# fasting insulin (up to 47,806 men/50,404 women),
file.out <- file.path(DIR_OUT, "FI_female_Lagou2018.txt")
N_study <- 50404 # FEMALE SAMPLE SIZE
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/FI_FG_Lagou2018/FI_STAGE1_2_3_SEX_GWAS_2018.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=snp, P=female_pvalue)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

# ======================================================================= #
# ============================  T2D_70kforT2D_Guarch2018 (SNPsnap mapped for rsIDs, ~8.5M SNPs)  =============================== #
# ======================================================================= #

# chrposID        chr_orig        position        Effect_Allele   Non_Effect_Allele       Freq_Effect_Allele      MAF     N_Effective_SampleSize  Effect  StdErr  Pvalue  Direction_Effects_Cohorts       Het_ISq Het_Pvalue      rsID    snp_maf chr     pos
# chr1:206676331  1       206676331       a       g       0.4827  0.4827  66939.7 -0.0301 0.0158  0.05655 +-----  0.0     0.9042  rs3860295       0.492   1       206676331
# chr1:38541140   1       38541140        a       t       0.0617  0.0617  66941.0 0.0045  0.0326  0.8897  -+--++  0.0     0.5289  rs12134181      0.04672 1       38541140
# chr1:233869132  1       233869132       t       g       0.476   0.476   66940.9 -0.0065 0.0157  0.6771  -++--+  26.7    0.2345  rs593879        0.4811  1       233869132



file.out <- file.path(DIR_OUT, "T2D_70kforT2D_Guarch2018.txt")
N_study <- 12931+57196
# 12,931 cases and 57,196 controls from European descent populations.
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_70kforT2D_Guarch2018/summarystatistics_metal_1KGPhase1_UK10K_imputed_PublicT2DGWAS_data_HetISq_75.w_snpsnap_rsid.txt.gz"
# GWAS of full sample including 113,006 individuals (32,384 cases and 80,622 controls)
df <- read_tsv(file.in)
head(df)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% rename(P=Pvalue, SNP=rsID) # rename
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  SCZ_Pardinas2018 (~8.1M with rsIDs)  =============================== #
# ======================================================================= #

### SNIPPET: *OBS* notice that the SNPids are pasted together with chr,pos and alleles
# SNP     Freq.A1 CHR     BP      A1      A2      OR      SE      P
# 10:100968448:T:AA       0.3519  10      100968448       t       aa      1.0024  0.01    0.812
# 10:101574552:A:ATG      0.4493  10      101574552       a       atg     0.98906 0.0097  0.2585
# 10:10222597:AT:A        0       10      10222597        a       at      0.9997  0.01    0.9777
# 10:102244152:A:AG       0.2008  10      102244152       a       ag      0.996705        0.0114  0.7701
# 10:102290698:A:AG       0.1899  10      102290698       a       ag      0.993024        0.0114  0.5431
### # ----> also contains rsIDs. e.g.:
# rs10000909:157367980:G:T        0.0273  105318
# rs10000910:166076842:A:C        0.4887  105318
# rs10000911:144136193:T:G        0.6327  105318
# rs10000912:3669829:A:G  0.1428  105318
# rs10000914:90547195:G:A 0.2647  105318
# rs10000916:184208957:C:T        0.7122  105318

file.out <- file.path(DIR_OUT, "SCZ_Pardinas2018.txt")
N_study <- 105318
# CASES = 40,675
# CONTROLS = 64,643
# TOTAL = 105,318
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SCZ_Pardinas2018/clozuk_pgc2.meta.sumstats.txt.gz"
# GWAS of full sample including 113,006 individuals (32,384 cases and 80,622 controls)
df <- read_tsv(file.in)
head(df)
df.filter <- df %>% filter(grepl("rs", SNP)) # keep only snps with rsIDs | only removes a few SNPs
df.filter.clean <- df.filter
df.filter.clean$SNP <- stringr::str_split_fixed(df.filter$SNP, pattern=":", n=Inf)[,1]
# df.tmp <- df %>% head %>% filter(grepl("rs", SNP)) %>% mutate(SNP1=stringr::str_split_fixed(SNP, pattern=":", n=Inf)[,1]) # *OBS SPECIFIC FOR SCZ_Pardinas2018.txt* cleaning SNP column
# ^ DOES NOT WORK - don't know why: "Error in mutate_impl(.data, dots) : Evaluation error: subscript out of bounds."
head(df.filter.clean)
df <- df.filter.clean %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

# ======================================================================= #
# ============================  INSOMNIA_Hammerschlag2017  =============================== #
# ======================================================================= #


# SNP     UNIQUE_ID       CHR     POS     A1      A2      EAF     OR      BETA    SE      P       N_analyzed      N_cases N_controls      INFO
# rs4881595       10:39145390_A_G 10      39145390        G       A       0.511773        1.00065 -0.00170299     0.00955389      0.85853 113006  32384   80622   0.967944
# rs375093451     10:39145519_G_C 10      39145519        C       G       0.0588732       1.00336 0.00345489      0.0202144       0.864293        113006  32384   80622   0.970578
# .       10:39145559_TGAATGGAATGGACTC_T  10      39145559        T       TGAATGGAATGGACTC        0.478942        1.00178 -0.000391588    0.010096        0.969063        113006  32384   80622   0.867022
# rs566496929     10:39146077_ATGGAATGGAATCGAATGGAATGGAC_A        10      39146077        A       ATGGAATGGAATCGAATGGAATGGAC      0.0472679       1.00616 0.00646568      0.0240618       0.788152        113006  32384   80622   0.878034
# rs200485340     10:39146089_C_G 10      39146089        G       C       0.448019        0.999235        -0.00304358     0.00960452      0.751326        113006  32384   80622   0.96705
# rs562508083     10:39146099_G_T 10      39146099        T       G       0.0458005       1.00707 0.00749806      0.0242236       0.756914        113006  32384   80622   0.89336

file.out <- file.path(DIR_OUT, "INSOMNIA_Hammerschlag2017.txt")
N_study <- 32384+80622
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INSOMNIA_Hammerschlag2017/Hammerschlag_NatGenet2017_insomnia_sumstats-full_090617.txt.gz"
# GWAS of full sample including 113,006 individuals (32,384 cases and 80,622 controls)
df <- read_tsv(file.in)
head(df)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  BIP_SCZ_PGC2018  =============================== #
# ======================================================================= #


# CHR     SNP     BP      A1      A2      FRQ_A_20129     FRQ_U_21524     INFO    OR      SE      P       Direction       HetPVa
# 8       rs62513865      101592213       T       C       0.0751  0.0775  0.923   0.96938 0.029   0.2837  --+--++++-++------+++-----+--+++        0.4506
# 8       rs79643588      106973048       A       G       0.0909  0.0951  0.997   0.95495 0.0253  0.06893 ----+--+++-++--+--+-++++--+-----        0.946
# 8       rs17396518      108690829       T       G       0.562   0.566   0.962   0.97629 0.0151  0.1126  +--+-+--++--+-+-+-+-+++--++---+-        0.3695
# 8       rs983166        108681675       A       C       0.558   0.563   0.993   0.97190 0.0149  0.05502 +--+-+-+++----+---+---+--++-+-++        0.4169



#### BDvsCONT
file.out <- file.path(DIR_OUT, "BIP_SCZ_BDvsCONT_PGC2018.txt")
N_study <- 20129+21524
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BIP_SCZ_PGC2018/BDvsCONT.sumstats.gz"
# BDvsCONT.sumstats.gz: Bipolar disorder cases (n=20,129) against an independent bipolar specific set of controls (n=21,524)
df <- read_tsv(file.in)
head(df)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


### SCZvsBD
file.out <- file.path(DIR_OUT, "BIP_SCZ_SCZvsBD_PGC2018.txt")
N_study <- 23585+15270
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BIP_SCZ_PGC2018/SCZvsBD.sumstats.gz"
# SCZvsBD.sumstats.gz: A subset of schizophrenia cases (n=23,585) and bipolar disorder cases (n=15,270) matched for ancestry and genotyping array platform.
df <- read_tsv(file.in)
head(df)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)



# ======================================================================= #
# ============================  AN_PGC_Duncan2017 (contains 10M SNPs, but only ~1M with rsIDs)  =============================== #
# ======================================================================= #


file.out <- file.path(DIR_OUT, "AN_PGC_Duncan2017.txt")
N_study <- 3495+10982 # --> 14477
# N_CASES = 3,495
# N_CONTROL = 10,982
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AN_PGC_Duncan2017/pgc.ed.freeze1.summarystatistics.July2017.txt.gz"
df <- read_tsv(file.in)
head(df)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  ADHD_PGC_Demontis2017  =============================== #
# ======================================================================= #

# CHR     SNP     BP      A1      A2      INFO    OR      SE      P
# 1       rs202152658     751343  A       T       0.884   1.03118 0.0221  0.1654
# 1       rs143225517     751756  T       C       0.884   0.96977 0.0221  0.1654
# 1       rs3094315       752566  A       G       0.943   0.96002 0.019   0.03208
# 1       rs3131972       752721  A       G       0.939   1.04081 0.0188  0.03301
# 1       rs3131971       752894  T       C       0.94    1.04300 0.0191  0.02772

file.out <- file.path(DIR_OUT, "ADHD_PGC_Demontis2017.txt")
N_study <- 19099+34194 # 53293
# adhd_eur_jun2017.gz: European ancestry meta-analysis (19,099 cases, 34,194 controls)
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/ADHD_PGC_Demontis2017/adhd_eur_jun2017.gz"
df <- read_tsv(file.in)
head(df)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  ASD_iPSYCH_PGC_Grove2018  =============================== #
# ======================================================================= #


# CHR     SNP     BP      A1      A2      INFO    OR      SE      P
# 8       rs62513865      101592213       T       C       0.949   1.00652 0.027   0.8086
# 8       rs79643588      106973048       A       G       0.997   1.01786 0.024   0.4606
# 8       rs17396518      108690829       T       G       0.987   0.96127 0.014   0.004647

file.out <- file.path(DIR_OUT, "ASD_iPSYCH_PGC_Grove2018.txt")
N_study <- 46350 # MAGMA docs: the total sample size, also when using case-control analysis results
# N_CASES = 18381
# N_CONTROL = 27969
# TOTAL = 18381+27969=46350
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/ASD_iPSYCH_PGC_Grove2018/iPSYCH-PGC_ASD_Nov2017.txt.gz"
df <- read_tsv(file.in)
head(df)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  ASD_iPSYCH_PGC_Anney2017  =============================== #
# ======================================================================= #

# chr     bp_hg19 snp     a1      a2      or      lb95    ub95    effect  se      p       frq_a1  info    N       direction
# 1       729679  rs4951859       C       G       1.05    0.97    1.14    .0504   .0422   .2323   .17     .621    13574   --+-+++-++++-+
# 1       731718  rs142557973     T       C       0.98    0.89    1.07    -.0239  .0475   .6143   .88     .657    13574   ++-+---+----+-
# 1       734349  rs141242758     T       C       0.96    0.88    1.06    -.0357  .0479   .4556   .88     .661    13574   ++-+---+----+-
# 1       736289  rs79010578      A       T       1.02    0.92    1.12    .0164   .0482   .7345   .13     .613    13574   -++-++-++-++-+
# 1       751756  rs143225517     T       C       0.97    0.89    1.05    -.0336  .0407   .4093   .86     .8      13574   ++-+----+---+-
# 1       752566  rs3094315       A       G       0.99    0.92    1.06    -.0111  .0375   .7663   .83     .805    13574   ++-+----+---+-


file.out <- file.path(DIR_OUT, "ASD_iPSYCH_PGC_Anney2017.txt")
N_study <- 13574 # MAGMA docs: the total sample size, also when using case-control analysis results
# N_CASES = 6197
# N_CONTROL = 7377
# TOTAL = 6197+7377=13574
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/ASD_iPSYCH_PGC_Anney2017/daner_AUT_meta14_CEU_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.with_header.tsv.gz"
df <- read_tsv(file.in)
head(df)
df <- df %>% rename(SNP=snp, P=p)
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

# ======================================================================= #
# ============================  MDD  =============================== #
# ======================================================================= #

# CHR     SNP     BP      A1      A2      FRQ_A_59851     FRQ_U_113154    INFO    OR      SE      P       ngt     Direction       HetISqt HetDf   HetPVa  Nca     Nco     Neff
# 8       rs62513865      101592213       T       C       0.0747  0.0733  0.957   1.01461 0.0153  0.3438  0       ---+++  0.0     5       0.7906  59851   113154  69115.85
# 8       rs79643588      106973048       A       G       0.092   0.092   0.999   1.02122 0.0136  0.1231  0       ++-+++  0.0     5       0.6847  59851   113154  69115.85

file.out <- file.path(DIR_OUT, "MDD_PGC_Wray2018.txt")
N_study <- 153780 # MAGMA docs: the total sample size, also when using case-control analysis results
# N_CASES = 130664-75607=55057
# N_CONTROL = 330470-231747=98723
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/MDD_PGC_Wray2018/MDD2018_ex23andMe.gz"
df <- read_tsv(file.in)
df.problems <- problems(df)
df <- df %>% filter(complete.cases(.)) # OBS: remove NA values
### dim
# 13554550 # df obs. incl. 122 parsing failures
# 13554551 # original file (wc -l) (incl header)
# 13554489 # df.red
df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  EA2  =============================== #
# ======================================================================= #

# MarkerName      CHR     POS     A1      A2      EAF     Beta    SE      Pval
# rs11130222      3       49901060        A       T       0.5765  0.026   0.003   4.581e-25

file.out <- file.path(DIR_OUT, "EA2_Okbay2016.txt")
N_study <- 217568 # MAGMA docs: the total sample size, also when using case-control analysis results

file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA2_Okbay2016/EduYears_Main.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=MarkerName, P=Pval) %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  EA3  =============================== #
# ======================================================================= #

# MarkerName      CHR     POS     A1      A2      EAF     Beta    SE      Pval
# rs13090388      3       49391082        T       C       0.6905  0.027   0.002   4.32e-54
# rs7630869       3       49522543        T       C       0.6922  0.027   0.002   4.62e-54

file.out <- file.path(DIR_OUT, "EA3_Lee2018.txt")
N_study <- 766345 # MAGMA docs: the total sample size, also when using case-control analysis results

file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA3_Lee2018/GWAS_EduYears_excl23andMe.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=MarkerName, P=Pval) %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df %>% write_tsv(file.out)



# ======================================================================= #
# ============================  SCZ  =============================== #
# ======================================================================= #

file.out <- file.path(DIR_OUT, "SCZ_Ripke2014.txt")
N_study <- 150064 # MAGMA docs: the total sample size, also when using case-control analysis results

file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SCZ_Ripke2014/ckqny.scz2snpres.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=snpid, P=p) %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df %>% write_tsv(file.out)

# ======================================================================= #
# ============================  HEIGHT  =============================== #
# ======================================================================= #

# MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  b       SE      p       N
# rs4747841       A       G       0.551   -0.0011 0.0029  0.70    253213

file.out <- file.path(DIR_OUT, "HEIGHT_Wood2014.txt")

file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/HEIGHT_Wood2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=MarkerName, P=p)
df <- df %>% select(SNP, P, N) # filter
df %>% write_tsv(file.out)

# ======================================================================= #
# ============================  LDL  =============================== #
# ======================================================================= #

# SNP_hg18        SNP_hg19        rsid    A1      A2      beta    se      N       P-value Freq.A1.1000G.EUR
# chr10:10000135  chr10:9960129   rs4747841       a       g       0.0037  0.0052  89138.00        0.7158  0.4908
# chr10:10000265  chr10:9960259   rs4749917       c       t       0.0033  0.0052  89138.00        0.7748  0.4908

### LDL
file.out <- file.path(DIR_OUT, "LIPIDS_LDL_Willer2013.txt")
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_LDL.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=rsid, P='P-value')
df <- df %>% select(SNP, P, N) # filter
df <- df %>% mutate(N=as.integer(N)) # make sure N is integer
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

### HDL
file.out <- file.path(DIR_OUT, "LIPIDS_HDL_Willer2013.txt")
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_HDL.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=rsid, P='P-value')
df <- df %>% select(SNP, P, N) # filter
df <- df %>% mutate(N=as.integer(N)) # make sure N is integer
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


### TC
file.out <- file.path(DIR_OUT, "LIPIDS_TC_Willer2013.txt")
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/LIPIDS_Willer2013/jointGwasMc_TC.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=rsid, P='P-value')
df <- df %>% select(SNP, P, N) # filter
df <- df %>% mutate(N=as.integer(N)) # make sure N is integer
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  WHR  =============================== #
# ======================================================================= #




# MarkerName      Allele1 Allele2 FreqAllele1HapMapCEU    b       se      p       N
# rs10011200      C       G       0.5333  0.017   0.0043  0.00011 142475
# rs8051831       T       C       0.0917  0.034   0.0089  0.00011 138860


### WHR
file.out <- file.path(DIR_OUT, "WHR_Shungin2015.txt")
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHR_COMBINED_EUR.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=MarkerName, P='p')
df <- df %>% select(SNP, P, N) # filter
df <- df %>% mutate(N=as.integer(N)) # make sure N is integer
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)


### WHR_adjBMI [**SPECIAL**]
file.out <- file.path(DIR_OUT, "WHR_adjBMI_Shungin2015.txt")
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz"
df <- read_tsv(file.in, col_types=cols('N'=col_double())) # *OBS specific to GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz*
# problems(df)
df <- df %>% rename(SNP=MarkerName, P='p')
df <- df %>% select(SNP, P, N) # filter
df <- df %>% mutate(N=as.integer(N)) # make sure N is integer
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)




# ======================================================================= #
# ============================  RA  =============================== #
# ======================================================================= #

file.out <- file.path(DIR_OUT, "RA_Okada2014.txt")
N_study <- 58284 # MAGMA docs: the total sample size, also when using case-control analysis results
# 14,361 RA cases and 43,923 conrols --> 58284

# SNPID   Chr     Position(hg19)  A1      A2      OR(A1)  OR_95%CIlow     OR_95%CIup      P-val
# chr1:751343     1       751343  A       T       0.85    0.75    0.96    0.01
# chr1:751756     1       751756  T       C       1.17    1.04    1.33    0.01
# rs3094315       1       752566  A       G       1.14    1.03    1.26    0.0093

file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/RA_Okada2014/RA_GWASmeta_European_v2.txt.gz"
df <- read_tsv(file.in)
df <- df %>% rename(SNP=SNPID, P='P-val') %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df %>% write_tsv(file.out)


# ======================================================================= #
# ============================  YENGO2018 (BMI+HEIGHT) =============================== #
# ======================================================================= #


### HEIGHT_Yengo2018
file.out <- file.path(DIR_OUT, "HEIGHT_Yengo2018.txt")
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz"
df <- read_tsv(file.in)
problems(df) # n_problems = 58 ---> we ignore them (see below)
#df <- df %>% rename(SNP=rsID, P=pval)
#df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df <- df %>% filter(complete.cases(.)) # OBS: remove NA values
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

### BMI_Yengo2018
file.out <- file.path(DIR_OUT, "BMI_Yengo2018.txt")
file.in <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz"
df <- read_tsv(file.in)
problems(df) # n_problems = 1 ---> we ignore them (see below)
#df <- df %>% rename(SNP=rsID, P=pval)
#df <- df %>% mutate(N=N_study) # modify
df <- df %>% select(SNP, P, N) # filter
df <- df %>% filter(complete.cases(.)) # OBS: remove NA values
stopifnot(!any(is.na(df))) # stop if any NA
df %>% write_tsv(file.out)

# df.hla <- df %>% filter((CHR==6) & (POS > 25e6) & (POS < 34e6) ) # ---> N=4877

### SNIPPET
# CHR     POS     SNP     Tested_Allele   Other_Allele    Freq_Tested_Allele_in_HRS       BETA    SE      P       N
# 7       92383888        rs10    A       C       0.06431  0.0006 0.0044   8.9e-01        598895
# 12      126890980       rs1000000       A       G       0.2219   0.0001 0.0023   9.5e-01        689928
# 4       21618674        rs10000010      T       C       0.5086   0.0008 0.0020   6.7e-01        785319
# 4       1357325 rs10000012      C       G       0.8634   0.0034 0.0028   2.3e-01        692463
# 4       37225069        rs10000013      A       C       0.7708  -0.0051 0.0024   3.4e-02        687856
# 4       84778125        rs10000017      T       C       0.2284   0.0042 0.0024   7.3e-02        686123

### EXAMPLE ISSUE WITH "N" column [Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz , line number 447476]:
### Solution ---> we just throw it away, because only 1 SNP has this problem (for BMI).
### See problems(df)
# CHR     POS     SNP     Tested_Allele   Other_Allele    Freq_Tested_Allele_in_HRS       BETA    SE      P       N
# ...
# 447476 6       87686337        rs12210954      A       C       0.97381  0.0180 0.0065   5.3e-03        667867
# 447477 6       6176144 rs12210959      T       C       0.7384  -0.0014 0.0022   5.1e-01        6e+05 *OBS*-----> N is converted to float
# 447478 6       76545684        rs12210963      T       C       0.8489   0.0008 0.0027   7.6e-01        692546
