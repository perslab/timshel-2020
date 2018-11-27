#!/usr/bin/env python2.7

import os
import sys
import numpy as np
import pandas as pd
import time

import glob

import multiprocessing

import pdb

###################################### USAGE ######################################

# time python gwas_to_rolypoly_fmt-map_chrpos.py |& tee log.alkes_group_UKBB.gwas_to_rolypoly_fmt-map_chrpos.txt

###################################### FUNCTIONS ######################################

### Read chrpos mapping file
def read_collection(file_collection):
    """Function that reads tab seperated gzip collection file"""
    print "START: reading CSV file PRIM..."
    start_time = time.time()
    f_tab = open(file_collection, 'r')
    df_collection = pd.read_csv(f_tab, index_col=False, header=0, delimiter="\t", compression="gzip")
    f_tab.close()
    elapsed_time = time.time() - start_time
    print "END: read CSV file PRIM into DataFrame in %s s (%s min)" % (elapsed_time, elapsed_time/60)

    # df_collection.head()
    #     rsID    snp_maf chr pos
    # 0   10:10001753 0.07455 10  10001753
    # 1   10:10001794 0.41050 10  10001794
    # 2   10:100023489    0.10540 10  100023489
    # 3   10:100025128    0.45430 10  100025128
    # 4   10:10002975 0.01193 10  10002975
    return df_collection

    

def read_and_process_gwas_file(file_gwas, GWAS_COL_SPECS, file_out_prefix, log_transform_beta):
    
    ### Allowing for missing columns
    flag_missing_se = False
    flag_missing_pval = False
    if GWAS_COL_SPECS[COL_OUT_se] is None:
        flag_missing_se = True
        print "*** OBS *** : no 'se' column given. Output column will be set np.nan"
    if GWAS_COL_SPECS[COL_OUT_pval] is None:
        flag_missing_pval = True
        print "*** OBS *** : no 'pval' column given. Output column will be set np.nan"
        
    ### Filtering dict by removing key-value pairs with the value None.
    GWAS_COL_SPECS_in_file = dict((k,v) for k,v in GWAS_COL_SPECS.iteritems() if v is not None) # REF https://stackoverflow.com/a/2544761/6639640
    print GWAS_COL_SPECS
    print GWAS_COL_SPECS_in_file


    df_gwas = pd.read_csv(file_gwas, usecols=GWAS_COL_SPECS_in_file.values(), index_col=False, header=0, delimiter="\s+") # read
    #df_gwas = pd.read_csv(file_gwas, usecols=GWAS_COL_SPECS_in_file.values(), index_col=False, header=0, delimiter="\t") # read
    assert( len(GWAS_COL_SPECS_in_file.values()) == len(GWAS_COL_SPECS_in_file.values()) ) # checking assumption (see below)
    df_gwas = df_gwas.rename(columns={v: k for k, v in GWAS_COL_SPECS_in_file.iteritems()}) # rename cols. REF: https://pandas.pydata.org/pandas-docs/stable/basics.html#basics-rename
        # ^ Reversing dict for renaming. ASSUMPTION: that all VALUES should be unique. REF: https://stackoverflow.com/a/483833/6639640
    if flag_missing_se: # add column of NaNs
        df_gwas[COL_OUT_se] = np.nan
    if flag_missing_pval: # add column of NaNs
        df_gwas[COL_OUT_pval] = np.nan

    ### check if all betas are positive. If they are, it is a odds ratio (OR).
    flag_beta_is_odds_ratio = (df_gwas.loc[:,"beta"] > 0).all() # OBS: an OR should always be > 0 and not >= 0, but we could use >= 0 to accomodate LDSC Z-score input.

    if log_transform_beta: # do log transformation
        print "************ DOING LOG OF BETA COLUMN TRANSFORMATION *************"
        df_gwas.loc[:,"beta"] = np.log(df_gwas.loc[:,"beta"]) # we know that the column is named "beta" now
        # ^^ # np.log(0) returns -inf and RuntimeWarning: divide by zero encountered in log
    elif flag_beta_is_odds_ratio:
        print "************ AUTOMATIC DETECTION - DOING LOG OF BETA COLUMN TRANSFORMATION *************"
        df_gwas.loc[:,"beta"] = np.log(df_gwas.loc[:,"beta"]) # we know that the column is named "beta" now
        # ^^ np.log(0) returns -inf and RuntimeWarning: divide by zero encountered in log
        file_out_tmp = "log.beta_log_transformed.{}.txt".format(file_out_prefix)
        with open(file_out_tmp, 'w') as fh_out:
            fh_out.write("{}: beta was log transformed (case control trait?).\n".format(file_out_prefix))

    ### We replace np.inf with np.nan, so we can throw them out 
    ### e.g. if we have taking the log of the betas, and some beta is zero: np.log(0) returns -inf.
    df_gwas.replace([np.inf, -np.inf], np.nan, inplace=True) 
    
    
    ### Check if there is any NaN values
    n_nan_beta_values = df_gwas.loc[:,"beta"].isnull().sum()
    if n_nan_beta_values > 0:
        tmp_string1 = "*WARNING* {}: beta column contains N={} NaN values.".format(file_out_prefix, n_nan_beta_values)
        tmp_string2 = "Will remove any SNPs with NaN values. Dim before removing NaN values: {}".format(df_gwas.shape)
        df_gwas.dropna(axis='index', how='any', subset=["beta"], inplace=True) # drop rows where  values in 'beta' is nan.
        tmp_string3 = "Done. Dim after removing NaN values: {}".format(df_gwas.shape)
        file_out_tmp = "log.nan_values_from_beta_col_removed.{}.txt".format(file_out_prefix)
        print tmp_string1
        print tmp_string2
        print tmp_string3
        with open(file_out_tmp, 'w') as fh_out:
            fh_out.write(tmp_string1+"\n")
            fh_out.write(tmp_string2+"\n")
            fh_out.write(tmp_string3+"\n")
    else:
        print "Wuhu, no NaN values found the beta column. All clear"
    
    return df_gwas


###################################### SCRIPT ######################################

def main(file_gwas, file_out_prefix):

    ### Read GWAS data
    print "START: reading GWAS..."
    df_gwas = read_and_process_gwas_file(file_gwas, GWAS_COL_SPECS, file_out_prefix, log_transform_beta=LOG_TRANSFORM_BETA)
    print "END: reading GWAS..."
    # df_gwas.head()
    #     rsID    beta    se  pval
    # 0   rs1000000   0.0001  0.0043  0.98140
    # 1   rs10000010  -0.0022 0.0029  0.43840
    # 2   rs10000012  -0.0096 0.0053  0.07009
    # 3   rs10000013  -0.0096 0.0043  0.02558
    # 4   rs10000017  -0.0038 0.0045  0.39840


    n_snps_not_found_in_snpsnap = sum(~df_gwas["rsID"].isin(df_collection["rsID"]))
    print "Number of SNPs in GWAS data: {}".format(len(df_gwas))
    print "Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: {}".format(n_snps_not_found_in_snpsnap)
    print "Percent SNPs not found: {:.2f} %".format(n_snps_not_found_in_snpsnap/float(len(df_gwas))*100)


    ### Join data
    # pd.merge(left, right, how='inner', on=None, left_on=None, right_on=None, left_index=False, right_index=False, sort=False, suffixes=('_x', '_y'), copy=True, indicator=False, validate=None)
    print "Joining data frames..."
    df_join = pd.merge(df_gwas, df_collection, how='inner', on="rsID")
    #     rsID    beta    se  pval    snp_maf chr pos
    # 0   rs1000000   0.0001  0.0043  0.98140 0.2237  12  126890980
    # 1   rs10000010  -0.0022 0.0029  0.43840 0.4901  4   21618674
    # 2   rs10000012  -0.0096 0.0053  0.07009 0.1402  4   1357325
    # 3   rs10000013  -0.0096 0.0043  0.02558 0.2296  4   37225069
    # 4   rs10000017  -0.0038 0.0045  0.39840 0.2475  4   84778125

    print "EXAMPLE output"
    print df_join.head()

    print "Dimensions of output file: {}".format(df_join.shape)

    ## Export
    print "START: exporting file..."
    file_out = "{}.gwassumstats.rolypoly_fmt.tab.gz".format(file_out_prefix)
    df_join.to_csv(file_out, sep="\t", index=False, compression='gzip')
    print "END: exported file: {}".format(file_out)


def main_parallel(file_gwas):
    file_out_prefix = os.path.basename(file_gwas).split(".")[0] # e.g. disease_RESPIRATORY_ENT.sumstats.gz --> disease_RESPIRATORY_ENT
    print "RUNNING {}".format(file_out_prefix)
    main(file_gwas, file_out_prefix)

def main_parallel2(i):
    file_gwas = list_file_gwas[i]
    file_out_prefix = list_file_output_prefix[i]
    print "RUNNING {}".format(file_out_prefix)
    main(file_gwas, file_out_prefix)


########################################################################################
###################################### GWAS PARAMETERS ######################################
########################################################################################

################# NULL_GWAS ##################

### **NB**: you must use white-space delimter to read NULL FILES: 
### df_gwas = pd.read_csv(file_gwas, usecols=GWAS_COL_SPECS_in_file.values(), index_col=False, header=0, delimiter="\s+") # read

### *OBS*: ORDER MUST BE THE SAME
list_file_gwas = ["/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/NULL_GWAS/plink_out.1KG_phase3_EUR_null_gwas.P{}.qassoc".format(i) for i in range(1,10+1)]

list_file_output_prefix = ["NULL_GWAS_1KG_phase3_EUR_N{}".format(i) for i in range(1,10+1)]

###  CHR                                  SNP         BP    NMISS       BETA         SE         R2        T            P
###    1                              1:11008      11008      503   -0.06232     0.1108  0.0006309  -0.5624       0.5741
###    1                              1:11012      11012      503   -0.06232     0.1108  0.0006309  -0.5624       0.5741
###    1                              1:13110      13110      503   -0.06193     0.1422  0.0003782  -0.4354       0.6635
###    1                          rs201725126      13116      503    0.05267    0.08318  0.0007997   0.6332       0.5269

COL_rsID = "SNP"
COL_beta = "BETA"
COL_se = "SE"
COL_pval = "P"
LOG_TRANSFORM_BETA=False #


################## body_HEIGHT_Wood2014 ##################

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/HEIGHT_Wood2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"
# file_out_prefix = "body_HEIGHT_Wood2014"

# ### MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  b       SE      p       N
# ### rs4747841       A       G       0.551   -0.0011 0.0029  0.70    253213
# ### rs4749917       T       C       0.436   0.0011  0.0029  0.70    253213

# COL_rsID = "MarkerName"
# COL_beta = "b"
# COL_se = "SE"
# COL_pval = "p"
# LOG_TRANSFORM_BETA=False # OBS --> some values are alrealy negative?



################## mental_AD_Lambert2013 ##################

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AD_Lambert2013/IGAP_stage_1.txt"
# file_out_prefix = "mental_AD_Lambert2013"

# ### Chromosome      Position        MarkerName      Effect_allele   Non_Effect_allele       Beta    SE      Pvalue
# ### 1       751343  rs28544273      A       T       -0.0146 0.0338  0.6651
# ### 1       751756  rs143225517     C       T       -0.0146 0.0338  0.6651
# ### 1       752566  rs3094315       G       A       -0.0122 0.0294  0.6773
# ### 1       753405  rs61770173      C       A       -0.0126 0.0339  0.7104

# COL_rsID = "MarkerName"
# COL_beta = "Beta"
# COL_se = "SE"
# COL_pval = "Pvalue"
# LOG_TRANSFORM_BETA=False # OBS --> some values are alrealy negative?


################## disease_RA_Okada2014 ##################

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/RA_Okada2014/RA_GWASmeta_European_v2.txt.gz"
# file_out_prefix = "disease_RA_Okada2014"

# ### SNPID   Chr     Position(hg19)  A1      A2      OR(A1)  OR_95%CIlow     OR_95%CIup      P-val
# ### chr1:751343     1       751343  A       T       0.85    0.75    0.96    0.01
# ### chr1:751756     1       751756  T       C       1.17    1.04    1.33    0.01
# ### rs3094315       1       752566  A       G       1.14    1.03    1.26    0.0093
# ### rs3131972       1       752721  A       G       0.88    0.79    0.97    0.009

# COL_rsID = "SNPID"
# COL_beta = "OR(A1)"
# COL_se = None
# COL_pval = "P-val"
# LOG_TRANSFORM_BETA=True # OBS


################## cov_INSOMNIA_Hammerschlag2017 ##################

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INSOMNIA_Hammerschlag2017/Hammerschlag_NatGenet2017_insomnia_sumstats-full_090617.txt.gz"
# file_out_prefix = "cov_INSOMNIA_Hammerschlag2017"

# ### SNP     UNIQUE_ID       CHR     POS     A1      A2      EAF     OR      BETA    SE      P       N_analyzed      N_cases N_controls      INFO
# ### rs4881595       10:39145390_A_G 10      39145390        G       A       0.511773        1.00065 -0.00170299     0.00955389      0.85853 113006  32384   80622   0.967944
# ### rs375093451     10:39145519_G_C 10      39145519        C       G       0.0588732       1.00336 0.00345489      0.0202144       0.864293        113006  32384   80622   0.970578
# ### .       10:39145559_TGAATGGAATGGACTC_T  10      39145559        T       TGAATGGAATGGACTC        0.478942        1.00178 -0.000391588    0.010096        0.969063        113006  32384   80622   0.867022
# ### rs566496929     10:39146077_ATGGAATGGAATCGAATGGAATGGAC_A        10      39146077        A       ATGGAATGGAATCGAATGGAATGGAC      0.0472679       1.00616 0.00646568      0.0240618       0.788152        113006  32384   80622   0.878034
# ### rs200485340     10:39146089_C_G 10      39146089        G       C       0.448019        0.999235        -0.00304358     0.00960452      0.751326        113006  32384   80622   0.96705
# ### rs562508083     10:39146099_G_T 10      39146099        T       G       0.0458005       1.00707 0.00749806      0.0242236       0.756914        113006  32384   80622   0.89336

# COL_rsID = "SNP"
# COL_beta = "BETA"
# COL_se = "SE"
# COL_pval = "P"
# LOG_TRANSFORM_BETA=False



################## cov_INTELLIGENCE_Sniekers2017 ##################

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/INTELLIGENCE_Sniekers2017/sumstats.txt.gz"
# file_out_prefix = "cov_INTELLIGENCE_Sniekers2017"

# ### Chromosome      position        rsid    ref     alt     MAF     Beta    SE      Zscore  p_value direction
# ### 1       100000012       rs10875231      T       G       0.234588        0.000298163293384453    0.00596326586768906     0.05    0.9599  +-++--+-
# ### 1       100000135       rs114947036     A       T       0.00180123      0.0410026320136857      0.0716829231008491      0.572   0.5676  ?-+?????
# ### 1       100000827       rs6678176       T       C       0.283746        0.00107059476923458     0.00560520821588785     0.191   0.8484  +-++----
# ### 1       100000843       rs78286437      T       C       0.0487127       -0.0104204709382541     0.0141198793201275      -0.738  0.4603  ?--?????
# ### 1       100000989       rs146963890     A       ATC     0.0465973       -0.00809008628857249    0.014420831173926       -0.561  0.575   ?--?????

# COL_rsID = "rsid"
# COL_beta = "Beta"
# COL_se = "SE"
# COL_pval = "p_value"
# LOG_TRANSFORM_BETA=False



################## BIP_PGC_2012 + MDD_PGC2013 ##################

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BIP_PGC2012/pgc.bip.full.2012-04.txt"
# file_out_prefix = "mental_BIP_PGC2012"

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/MDD_PGC2013/pgc.mdd.full.2012-04.txt"
# file_out_prefix = "mental_MDD_PGC2013"


### snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
### rs12565286      1       711153  C       G       0.9887  0.0821  0.8899  0.3693  0       0.0678
### rs11804171      1       713682  A       T       0.9887  0.0821  0.8899  0.3693  0       .
### rs2977670       1       713754  C       G       1.011   0.0822  0.8936  0.3684  0       .
### rs12138618      1       740098  A       G       1.1241  0.0927  0.2069  0.2765  0       0.05833
### rs3094315       1       742429  A       G       1.0037  0.0353  0.9175  0.7758  0       0.8448
### rs3131968       1       744055  A       G       1.0021  0.0416  0.9605  0.6888  0       0.125
### rs12562034      1       758311  A       G       1.0716  0.0579  0.2321  0.4646  0       0.09167

# COL_rsID = "snpid"
# COL_beta = "or"
# COL_se = "se"
# COL_pval = "pval"
# LOG_TRANSFORM_BETA=True # *OBS*



################## ASD_iPSYCH_PGC_Grove2018 ##################

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/ASD_iPSYCH_PGC_Grove2018/daner_AUT_meta14_CEU_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.with_header.tsv.gz"
# file_out_prefix = "mental_ASD_iPSYCH_PGC_Grove2018"

# ### chr     bp_hg19 snp     a1      a2      or      lb95    effect  se      p       frq_a1  info    direction
# ### 1       729679  rs4951859       C       G       1.05    0.97    1.14    .0504   .0422   .2323   .17     .621    13574   --+-+++-++++-+
# ### 1       731718  rs142557973     T       C       0.98    0.89    1.07    -.0239  .0475   .6143   .88     .657    13574   ++-+---+----+-
# ### 1       734349  rs141242758     T       C       0.96    0.88    1.06    -.0357  .0479   .4556   .88     .661    13574   ++-+---+----+-
# ### 1       736289  rs79010578      A       T       1.02    0.92    1.12    .0164   .0482   .7345   .13     .613    13574   -++-++-++-++-+

# ### chr ........... chromosome
# ### bp_hg19 ....... location on chromosome (based on reference hg19 +1)
# ### snp ........... variant identifier
# ### a1 ............ reference allele
# ### a2 ............ alternate allele
# ### or ............ odds ratio - derived from "effect", based on reference allele
# ### lb95 / ub95 ... lower and upper boundries (95%) CI around the reported odds ratio - derived from "se", based on reference allele
# ### effect ........ effect estimate from regression, based on reference allele
# ### se ............ standard error of the effect, based on reference allele
# ### p ............. p-value
# ### frq_a1 ........ approximate (rounded) reference (a1) allele frequency for all participants in the data set *not in noswm3 - created using bespoke pipeline*
# ### info .......... INFO score from imputation, ratio of variances, can exceed 1                               *not in noswm3 - created using bespoke pipeline*
# ### direction ..... the direction of effect in the contributing data sets (meta-analysis) (note// ? == missing)

# COL_rsID = "snp"
# COL_beta = "effect" # *OBS*: notice that there is also a "or" column. *I THINK 'or' could just as well be used*
# COL_se = "se"
# COL_pval = "p"
# LOG_TRANSFORM_BETA=True # *OBS*




################## SWB_Neuroticism_DS_Okbay2016 ##################

# ### *OBS*: ORDER MUST BE THE SAME
# list_file_gwas = ["/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/DS_Full.txt.gz",
# "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/Neuroticism_Full.txt.gz",
# "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/SWB_Full.txt.gz"]

# list_file_output_prefix = ["mental_DS_Okbay2016",
# "mental_Neuroticism_Okbay2016",
# "mental_SWB_Okbay2016"]


# # file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/DS_Full.txt.gz"
# # file_out_prefix = "mental_DS_Okbay2016"
# # file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/Neuroticism_Full.txt.gz"
# # file_out_prefix = "mental_Neuroticism_Okbay2016"
# # file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SWB_Neuroticism_DS_Okbay2016/SWB_Full.txt.gz"
# # file_out_prefix = "mental_SWB_Okbay2016"

# ### MarkerName      CHR     POS     A1      A2      EAF     Beta    SE      Pval
# ### rs2572431       8       11105077        T       C       0.5951  0.028   0.003   4.197e-16
# ### rs9286062       8       11118527        A       G       0.5784  0.027   0.003   2.32e-15
# ### rs2293855       8       11177410        A       G       0.3769  -0.028  0.004   3.781e-15

# COL_rsID = "MarkerName"
# COL_beta = "Beta" # OBS: this is the trick to LDSC
# COL_se = "SE" # this column is not present. 
# COL_pval = "Pval"
# LOG_TRANSFORM_BETA=False #




################## LDSC gwas_sumstats_ldsc/alkesgroup/sumstats_formatted ##################


# list_files = glob.glob("/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/alkesgroup/sumstats_formatted/*.sumstats")
# print list_files

# ### *OBS*: "elif flag_beta_is_odds_ratio:" was disabled for these traits.

# ### SNP A1  A2  N   CHISQ   Z
# ### rs4970383   A   C   14267.000   2.010724    1.418
# ### rs4475691   T   C   14267.000   1.236544    1.112
# ### rs1806509   C   A   14267.000   1.0816  1.040
# ### rs7537756   G   A   14267.000   0.521284    0.722
# ### rs13302982  A   G   14267.000   2.408704    1.552

# COL_rsID = "SNP"
# COL_beta = "Z" # OBS: this is the trick to LDSC
# COL_se = None # this column is not present. 
# COL_pval = None
# LOG_TRANSFORM_BETA=False #


#### *WARNING* Type_1_Diabetes: beta column contains N=102 NaN values.
#### Will remove any SNPs with NaN values. Dim before removing NaN values: (915517, 4)
#### Done. Dim after removing NaN values: (915415, 4)
#### *WARNING* IBD: beta column contains N=152 NaN values.
#### Will remove any SNPs with NaN values. Dim before removing NaN values: (1078060, 4)
#### Done. Dim after removing NaN values: (1077908, 4)

################## mental_AN_Boraska2014 ##################

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AN_Boraska2014/gcan_meta.out.gz"
# file_out_prefix = "mental_AN_Boraska2014"

# ### chromosome      position        SNP     reference_allele        other_allele    eaf     OR      OR_se   OR_95L  OR_95U  z       p_sanger        _-log10_p-value q_statistic     q_p-value       i2      n_studies       n_samples       effects
# ### 1       100208721       rs11166383      G       A       -9      1.082253        0.272185        0.548771        2.134356        0.228132        0.819536        0.086432        0.332114        0.564417        0.000000        2       -9      -??????+???????
# ### 1       100212478       rs11166384      T       C       -9      1.035175        0.043051        0.950794        1.127043        0.796884        0.425523        0.371077        22.540501       0.068170        0.378896        15      -9      ++-+--+++----++
# ### 1       100214368       rs6686328       G       A       -9      0.874292        0.119020        0.641013        1.192467        -0.848381       0.396226        0.402057        5.708951        0.126662        0.474509        4       -9      +?-???+??-?????

# COL_rsID = "SNP"
# COL_beta = "OR"
# COL_se = "OR_se"
# COL_pval = "p_sanger"
# LOG_TRANSFORM_BETA=True # *OBS*


################## mental_AN_Duncan2017 *NOT COMPLETD* ##################

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/AN_PGC_Duncan2017/pgc.ed.freeze1.summarystatistics.July2017.txt.gz"
# file_out_prefix = "mental_AN_Duncan2017"

# ### CHR     SNP     BP      A1      A2      INFO    OR      SE      P       ngt
# ### 1       chr1:100012367  100012367       T       C       0.842   3.84164 1.5     0.3694  0
# ### 1       chr1:100020280  100020280       T       C       0.668   0.00290 4.48    0.1922  0
# ### 1       chr1:100023954  100023954       C       G       0.748   3.18930 1.443   0.4216  0
# ### 1       chr1:100035002  100035002       G       GT      0.755   3.48860 1.918   0.5147  0

# COL_rsID = "SNP"
# COL_beta = "BETA"
# COL_se = "SE"
# COL_pval = "P"
# LOG_TRANSFORM_BETA=False

################## BMI_HEIGHT_Yengo2018 ##################

# ### BMI
# # file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz"
# # file_out_prefix = "body_BMI_Yengo2018"

# ### HEIGHT
# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz"
# file_out_prefix = "body_height_Yengo2018"


# ### CHR     POS     SNP     Tested_Allele   Other_Allele    Freq_Tested_Allele_in_HRS       BETA    SE      P       N
# ### 7       92383888        rs10    A       C       0.06431  0.0006 0.0044   8.9e-01        598895
# ### 12      126890980       rs1000000       A       G       0.2219   0.0001 0.0023   9.5e-01        689928
# ### 4       21618674        rs10000010      T       C       0.5086   0.0008 0.0020   6.7e-01        785319
# ### 4       1357325 rs10000012      C       G       0.8634   0.0034 0.0028   2.3e-01        692463
# ### 4       37225069        rs10000013      A       C       0.7708  -0.0051 0.0024   3.4e-02        687856


# COL_rsID = "SNP"
# COL_beta = "BETA"
# COL_se = "SE"
# COL_pval = "P"
# LOG_TRANSFORM_BETA=False


################## MENOPAUSE_Day2015 ##################

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/MENOPAUSE_Day2015/Menopause_HapMap2_DayNG2015_18112015.txt.gz"
# file_out_prefix = "repro_menopause_Day2015"

# ## MarkerName      allele1 allele2 HapMap_eaf      effect  stderr  p
# ## rs10    a       c       0.03    -0.01   0.05    9.0e-01
# ## rs1000000       a       g       0.37    -0.01   0.02    5.8e-01
# ## rs10000010      t       c       0.57    0.01    0.02    6.7e-01
# ## rs10000012      c       g       0.81    -0.07   0.03    1.9e-02

# COL_rsID = "MarkerName"
# COL_beta = "effect"
# COL_se = "stderr"
# COL_pval = "p"
# LOG_TRANSFORM_BETA=False

# Plus_iCOGs_Beta* =  SNP effect size for the GWAS meta-analysis including iCOGs samples (max N = 182,416)
# Plus_iCOGs_P* =  SNP p-value for the GWAS meta-analysis including iCOGs samples
# *this will give same estimates as GWAS_Beta/P if SNP is not present in iCOGs.


################## MENARCHE_Perry2014 ##################

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/MENARCHE_Perry2014/Menarche_Nature2014_GWASMetaResults_17122014.txt.gz"
# file_out_prefix = "repro_menarche_Perry2014"

# ## MarkerName      Effect_Allele   Other_Allele    HapMap_Freq     GWAS_Beta       GWAS_P  iCOGs   Plus_iCOGs_Beta Plus_iCOGs_P
# ## rs10754651      t       c       .02     -.17    6.4e-04 N       -.17    6.4e-04
# ## rs7751829       t       c       .98     -.16    3.0e-07 Y       -.089   6.1e-05
# ## rs10428740      t       c       .01     -.15    1.0e-02 N       -.15    1.0e-02
# ## rs9895917       t       c       .02     -.15    7.8e-04 N       -.15    7.8e-04

# COL_rsID = "MarkerName"
# COL_beta = "Plus_iCOGs_Beta"
# COL_se = None # this column is not present. 
# COL_pval = "Plus_iCOGs_P"
# LOG_TRANSFORM_BETA=False

## Plus_iCOGs_Beta* =  SNP effect size for the GWAS meta-analysis including iCOGs samples (max N = 182,416)
## Plus_iCOGs_P* =  SNP p-value for the GWAS meta-analysis including iCOGs samples
## *this will give same estimates as GWAS_Beta/P if SNP is not present in iCOGs.

################## disease_T2D_Guarch2018 ##################

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/T2D_70kforT2D_Guarch2018/summarystatistics_metal_1KGPhase1_UK10K_imputed_PublicT2DGWAS_data_HetISq_75.w_snpsnap_rsid.txt.gz"
# file_out_prefix = "disease_T2D_Guarch2018"

# COL_rsID = "rsID"
# COL_beta = "Effect" # these also contain negative values, so they are likely already log-transformed
# COL_se = "StdErr"
# COL_pval = "Pvalue"
# LOG_TRANSFORM_BETA=False

################## EA3 Lee2018 ##################
# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA3_Lee2018/GWAS_EduYears_excl23andMe.txt.gz"
# file_out_prefix = "cov_EduYears_Lee2018"

# ## MarkerName      CHR     POS     A1      A2      EAF     Beta    SE      Pval
# ## rs13090388      3       49391082        T       C       0.6905  0.027   0.002   4.32e-54
# ## rs7630869       3       49522543        T       C       0.6922  0.027   0.002   4.62e-54
# ## rs7623659       3       49414791        T       C       0.3095  0.027   0.002   4.77e-54

# COL_rsID = "MarkerName"
# COL_beta = "Beta"
# COL_se = "SE"
# COL_pval = "Pval"
# LOG_TRANSFORM_BETA=False


################## EA2 Okbay2016 ##################
# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/EA2_Okbay2016/EduYears_Main.txt.gz"
# file_out_prefix = "cov_EduYears_Okbay2016"

# COL_rsID = "MarkerName"
# COL_beta = "Beta"
# COL_se = "SE"
# COL_pval = "Pval"
# LOG_TRANSFORM_BETA=False

################## WHR ##################
## file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHR_COMBINED_AllAncestries.txt.gz"
## file_out_prefix = "body_WHR_Shungin2015_ALL"

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz"
# file_out_prefix = "body_WHRadjBMI_Shungin2015_EUR"

# file_gwas = "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/WHR_Shungin2015/GIANT_2015_WHR_COMBINED_EUR.txt.gz"
# file_out_prefix = "body_WHR_Shungin2015_EUR"


# ### MarkerName      Chr     Pos     Allele1 Allele2 FreqAllele1HapMapCEU    b       se      p       N
# ### rs13248326      8       116842464       G       A       0.075   0.036   0.0094  0.0001  111586
# ### rs10842702      12      26331326        C       G       0.3     0.014   0.0036  0.0001  215116
# ### rs3843351       3       66466285        C       T       0.9083  0.033   0.0084  0.0001  132858
# ### rs8100151       19      38496239        C       G       0.8417  0.025   0.0065  0.0001  132505
# ### rs13274892      8       23641999        G       C       0.7417  0.019   0.005   0.0001  144587
# ### rs2289921       11      8453187 C       G       0.5917  0.016   0.0042  0.0001  144595


# COL_rsID = "MarkerName"
# COL_beta = "b"
# COL_se = "se"
# COL_pval = "p"
# LOG_TRANSFORM_BETA=False

################## Lipids ##################
# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/jointGwasMc_HDL.txt.gz"
# file_out_prefix = "lipids_HDL_Willer2013"
# Number of SNPs in GWAS data: 2437751
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 63826
# Percent SNPs not found: 2.62 %

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/jointGwasMc_LDL.txt.gz"
# file_out_prefix = "lipids_LDL_Willer2013"
# Number of SNPs in GWAS data: 2447441
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 66167
# Percent SNPs not found: 2.70 %

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/jointGwasMc_TC.txt.gz"
# file_out_prefix = "lipids_TC_Willer2013"

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/jointGwasMc_TG.txt.gz"
# file_out_prefix = "lipids_TG_Willer2013"

# mv jointGwasMc_HDL.gwassumstats.rolypoly_fmt.tab.gz lipids_HDL_Willer2013.gwassumstats.rolypoly_fmt.tab.gz
# mv jointGwasMc_LDL.gwassumstats.rolypoly_fmt.tab.gz lipids_LDL_Willer2013.gwassumstats.rolypoly_fmt.tab.gz
# mv jointGwasMc_TC.gwassumstats.rolypoly_fmt.tab.gz lipids_TC_Willer2013.gwassumstats.rolypoly_fmt.tab.gz
# mv jointGwasMc_TG.gwassumstats.rolypoly_fmt.tab.gz lipids_TG_Willer2013.gwassumstats.rolypoly_fmt.tab.gz

# list_files = glob.glob("/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/jointGwasMc_*.txt.gz")
# print list_files
#     # /Users/djw472/data/GWAS-sumstats/alkesgroup-collection/UKBB/body_BMIz.sumstats.gz
#     # /Users/djw472/data/GWAS-sumstats/alkesgroup-collection/UKBB/body_HEIGHTz.sumstats.gz
# # file_out_prefix = <NOT DEFINED FOR ALKES / LOOP MODE>

# COL_rsID = "rsid"
# COL_beta = "beta"
# COL_se = "se"
# COL_pval = "P-value"

################## SCZ_Ripke2014 ##################
# hg19chrc        snpid   a1      a2      bp      info    or      se      p       ngt
# chr1    rs4951859       C       G       729679  0.631   0.97853 0.0173  0.2083  0
# chr1    rs142557973     T       C       731718  0.665   1.01949 0.0198  0.3298  0
# chr1    rs141242758     T       C       734349  0.666   1.02071 0.02    0.3055  0

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/SCZ_Ripke2014/ckqny.scz2snpres.gz"
# file_out_prefix = "mental_SCZ_Ripke2014"

# COL_rsID = "snpid"
# COL_beta = "or"
# COL_se = "se"
# COL_pval = "p"
# LOG_TRANSFORM_BETA=True # do logtransformation


# Number of SNPs in GWAS data: 9444230
# Number of SNPs in GWAS data *NOT found* in SNP chr pos mapping file: 1788142
# Percent SNPs not found: 18.93 %
# Joining data frames...



################## alkes_lo_UKBB ##################

# SNP            CHR      POS     A1      A2      REF     EAF             Beta            se              P       N       INFO
# rs10399793      1       49298   T       C       T       0.37622         0.000334086     0.000731539     7.0E-01 459324  0.342797
# rs2462492       1       54676   C       T       C       0.599409        -0.000937692    0.000724671     1.7E-01 459324  0.340158
# rs3107975       1       55326   T       C       T       0.991552        0.00706826      0.0040343       8.3E-02 459324  0.324228
# rs74447903      1       57033   T       C       T       0.998221        0.00919789      0.00897642      3.8E-01 459324  0.296256
# 1:70728_C_T     1       70728   C       T       C       0.997834        0.00859618      0.00730365      2.6E-01 459324  0.365713
# rs2462495       1       79033   A       G       A       0.00129115      0.00321513      0.00929391      7.7E-01 459324  0.536566
# rs114608975     1       86028   T       C       T       0.896384        0.000549046     0.00115835      5.8E-01 459324  0.340885

# list_files = ["/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/cov_SMOKING_STATUS.sumstats.gz"] # SINGLE-TEST
# # list_files = glob.glob("/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/alkesgroup_UKBB/*sumstats.gz") # get all files
# print list_files

# ## file_out_prefix = <NOT DEFINED FOR ALKES / LOOP MODE>

# COL_rsID = "SNP"
# COL_beta = "Beta"
# COL_se = "se"
# COL_pval = "P"
# LOG_TRANSFORM_BETA=False # do logtransformation

################## BMI_Locke2015 ##################

# file_gwas = "/Users/djw472/data/GWAS-sumstats/timshel-collection/BMI_Locke2015/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq" # (cannot read .zip from GUI OSX zipped file)
# file_out_prefix = "body_BMI_Locke2015_All"

# file_gwas = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/BMI_Locke2015/EUR_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz" #
# file_out_prefix = "body_BMI_Locke2015_EUR"

# COL_rsID = "SNP"
# COL_beta = "b"
# COL_se = "se"
# COL_pval = "p"
# LOG_TRANSFORM_BETA=False # do logtransformation

#####################################################################################
###################################### DROPPED ######################################
#####################################################################################

################## BMI.Finucane_UKB ##################
# file_gwas = "/Users/djw472/data/GWAS-sumstats/finucane-collection/public_sumstats/PASS_BMI1.sumstats"
# file_out_prefix = "BMI.Finucane_UKB"
# # SNP     A1      A2      N       CHISQ   Z
# COL_rsID = "SNP"
# COL_beta = "Z"
# COL_se = "se"
# COL_pval = "p"


#####################################################################################
###################################### CONSTANTS ######################################
#####################################################################################


# file_collection = "/Users/djw472/Dropbox/0_Projects/p_sc_genetics/analysis/src/snpsnap_EUR_1KG_phase3-chrpos_mapping.tab.gz" # OSX
file_collection = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/snp_mapping/snpsnap_EUR_1KG_phase3-chrpos_mapping.tab.gz" # Ygg

COL_OUT_rsID = "rsID"  # do not change unless you really mean it - used to match rolypoly names
COL_OUT_beta = "beta" # do not change unless you really mean it - used to match rolypoly names
COL_OUT_se = "se" # do not change unless you really mean it - used to match rolypoly names
COL_OUT_pval = "pval" # do not change unless you really mean it - used to match rolypoly names
# GWA_OUTS_COL_SPECSCOL_rsIDCOL_rsID:COL_OUT_beta,COL_beta:COL_OUT_se,COL_se:COL_OUT_pval,COL_pval:"pval"} # do not change this. This "direction" of the dict is most convenient for pandas col rename
GWAS_COL_SPECS = {COL_OUT_rsID:COL_rsID,
                COL_OUT_beta:COL_beta,
                COL_OUT_se:COL_se,
                COL_OUT_pval:COL_pval} # do not change this. This "direction" of the dict makes it possible to work with multiple missing columns (None)

#####################################################################################
######################################## MAIN  ######################################
#####################################################################################

df_collection = read_collection(file_collection)

### Normal mode
# main(file_gwas, file_out_prefix)

### Alkes/loop mode
# for file_gwas in list_files:
#     file_out_prefix = os.path.basename(file_gwas).split(".")[0]
#     print "RUNNING {}".format(file_out_prefix)
#     main(file_gwas, file_out_prefix)

### Alkes/loop mode - parallel
### REMEMBER TO MAKE list_files 
# print "Starting pool..."
# pool = multiprocessing.Pool(35)
# pool.map(main_parallel, list_files)

## Loop mode - parallel2 | USE WEHN list_file_output_prefix is DEFINED
print "Starting pool..."
pool = multiprocessing.Pool(10)
pool.map(main_parallel2, range(len(list_file_gwas)))




print "SCRIPT ENDED"

