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


# body_BMI_Yengo2018      2322004 1305090 56.2053295343   (1016914, 7)
# body_BMIz       7594114 6415733 84.4829693102   (1178381, 7)
# body_HEIGHTz    7594110 6415729 84.482961137    (1178381, 7)
# body_WHRadjBMIz 7594107 6415726 84.4829550071   (1178381, 7)
# cov_EDU_COLLEGE 7594090 6415712 84.4829597753   (1178378, 7)
# cov_EDU_YEARS   7594090 6415712 84.4829597753   (1178378, 7)
# disease_T2D     7594109 6415728 84.4829590937   (1178381, 7)
# mental_NEUROTICISM      7593928 6415560 84.482760437    (1178368, 7)
# repro_MENARCHE_AGE      7594294 6415911 84.4833107594   (1178383, 7)
# lipids_TG_Willer2013    2375076 1342743 56.5347382568   (1032333, 7)
# cov_EduYears_Okbay2016  7719283 6548897 84.838151419    (1170386, 7)
# body_BMI_Locke2015_All  2466089 1390631 56.3901383932   (1075458, 7)
# cov_EduYears_Lee2018    7547186 6371203 84.4182586728   (1175983, 7)
# mental_SCZ_Ripke2014    7656088 6481934 84.6637865187   (1174154, 7)
# body_WHR_Shungin2015    2451735 1379916 56.2832443147   (1071819, 7)
# blood_EOSINOPHIL_COUNT  7594109 6415730 84.4829854299   (1178379, 7)
# blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT     7594077 6415698 84.4829200441   (1178379, 7)
# blood_LYMPHOCYTE_COUNT  7594102 6415722 84.4829579587   (1178380, 7)
# blood_MEAN_CORPUSCULAR_HEMOGLOBIN       7594097 6415717 84.4829477422   (1178380, 7)
# blood_MEAN_PLATELET_VOL 7594113 6415733 84.482980435    (1178380, 7)
# blood_MEAN_SPHERED_CELL_VOL     7594066 6415688 84.4829107358   (1178378, 7)
# blood_MONOCYTE_COUNT    7594119 6415739 84.4829926947   (1178380, 7)
# blood_PLATELET_COUNT    7594105 6415725 84.4829640886   (1178380, 7)
# blood_PLATELET_DISTRIB_WIDTH    7594106 6415726 84.4829661319   (1178380, 7)
# blood_RBC_DISTRIB_WIDTH 7594112 6415732 84.4829783917   (1178380, 7)
# blood_RED_COUNT 7594108 6415728 84.4829702185   (1178380, 7)
# blood_WHITE_COUNT       7594103 6415723 84.482960002    (1178380, 7)
# bmd_HEEL_TSCOREz        7594074 6415701 84.4829929232   (1178373, 7)
# body_BALDING1   7593699 6415347 84.4825031911   (1178352, 7)
# body_BALDING4   7593699 6415347 84.4825031911   (1178352, 7)
# bp_DIASTOLICadjMEDz     7594128 6415747 84.4829979163   (1178381, 7)
# bp_SYSTOLICadjMEDz      7594128 6415747 84.4829979163   (1178381, 7)
# cov_SMOKING_STATUS      7594108 6415728 84.4829702185   (1178380, 7)
# disease_AID_ALL 7594109 6415728 84.4829590937   (1178381, 7)
# disease_ALLERGY_ECZEMA_DIAGNOSED        7594102 6415721 84.4829447906   (1178381, 7)
# disease_ASTHMA_DIAGNOSED        7594102 6415721 84.4829447906   (1178381, 7)
# disease_CARDIOVASCULAR  7594109 6415728 84.4829590937   (1178381, 7)
# disease_DERMATOLOGY     7594109 6415728 84.4829590937   (1178381, 7)
# disease_HI_CHOL_SELF_REP        7594109 6415728 84.4829590937   (1178381, 7)
# disease_HYPOTHYROIDISM_SELF_REP 7594109 6415728 84.4829590937   (1178381, 7)
# disease_RESPIRATORY_ENT 7594109 6415728 84.4829590937   (1178381, 7)
# impedance_BASAL_METABOLIC_RATEz 7594109 6415728 84.4829590937   (1178381, 7)
# lung_FEV1FVCzSMOKE      7594072 6415689 84.4828571549   (1178383, 7)
# lung_FVCzSMOKE  7594072 6415689 84.4828571549   (1178383, 7)
# other_MORNINGPERSON     7594114 6415730 84.4829298059   (1178384, 7)
# pigment_HAIR_blonde     7594109 6415728 84.4829590937   (1178381, 7)
# pigment_HAIR_darkbrown  7594109 6415728 84.4829590937   (1178381, 7)
# pigment_HAIR    7594129 6415751 84.4830394638   (1178378, 7)
# pigment_SKIN    7594086 6415708 84.4829516021   (1178378, 7)
# pigment_SUNBURN 7594169 6415788 84.4830816907   (1178381, 7)
# pigment_TANNING 7594099 6415719 84.4829518288   (1178380, 7)
# repro_MENOPAUSE_AGE     7594067 6415693 84.4829654518   (1178374, 7)
# lipids_HDL_Willer2013   2381274 1346971 56.5651411807   (1034303, 7)
# lipids_LDL_Willer2013   2373925 1341910 56.5270596165   (1032015, 7)
# lipids_TC_Willer2013    2381020 1346785 56.5633636005   (1034235, 7)
# body_BMI_Locke2015_EUR  2465956 1390564 56.3904627658   (1075392, 7)
# disease_T2D_Guarch2018  8675737 7501111 86.4607928986   (1174626, 7)
# repro_menarche_Perry2014        2394306 1340603 55.9912976871   (1053703, 7)
# repro_menopause_Day2015 2378000 1330427 55.9473086627   (1047573, 7)
# body_height_Yengo2018   2307857 1297346 56.2143148384   (1010511, 7)
# mental_AN_Boraska2014   1135695 96970   8.53838398514   (1038725, 7)
# body_WHRadjBMI_Shungin2015_EUR  2433188 1362669 56.0034407534   (1070519, 7)
# body_WHR_Shungin2015_EUR        2451665 1379865 56.2827710964   (1071800, 7)
# cov_INSOMNIA_Hammerschlag2017   8738505 7566048 86.5828651468   (1172457, 7)
# cov_INTELLIGENCE_Sniekers2017   8735209 7563025 86.580927829    (1172184, 7)
# mental_AD_Lambert2013   6712749 5574216 83.0392436839   (1138533, 7)
# mental_ASD_iPSYCH_PGC_Grove2018 5639435 4567116 80.9853469364   (1072319, 7)
# mental_BIP_PGC2012      2381654 1327946 55.7573014384   (1053708, 7)
# mental_DS_Okbay2016     6405799 5285275 82.507662198    (1120524, 7)
# mental_MDD_PGC2013      1172105 102032  8.70502216098   (1070073, 7)
# mental_Neuroticism_Okbay2016    6405747 5285223 82.5075202002   (1120524, 7)
# mental_SWB_Okbay2016    2245164 1232819 54.909975396    (1012345, 7)
# body_HEIGHT_Wood2014    2459357 1386784 56.3880721668   (1072573, 7)
# NULL_GWAS_1KG_phase3_EUR_N10    9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N1     9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N2     9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N3     9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N4     9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N5     9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N6     9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N7     9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N8     9535059 8357442 87.6496097192   (1177617, 7)
# NULL_GWAS_1KG_phase3_EUR_N9     9535059 8357442 87.6496097192   (1177617, 7)


###################################### FUNCTIONS ######################################

### Read chrpos mapping file
def read_hapmap_snplist():
    file_hapmap = "/projects/timshel/sc-genetics/ldsc/data/w_hm3.snplist"

    ### SNIPPET
    # SNP     A1      A2
    # rs3094315       G       A
    # rs3131972       A       G
    # rs3131969       A       G
    # rs1048488       C       T
    # rs3115850       T       C


    print "START: reading file Hapmap..."
    start_time = time.time()
    f_tab = open(file_hapmap, 'r')
    df_collection = pd.read_csv(f_tab, index_col=False, usecols=['SNP'], header=0, delimiter="\t")
    f_tab.close()
    elapsed_time = time.time() - start_time
    print "END: read hapmap DataFrame in %s s (%s min)" % (elapsed_time, elapsed_time/60)

    df_collection = df_collection.rename(columns={"SNP": "rsID"}) # rename cols. REF: https://pandas.pydata.org/pandas-docs/stable/basics.html#basics-rename

    return df_collection



###################################### SCRIPT ######################################

def main(file_gwas, file_out_prefix):

    ### Read GWAS data
    print "START: reading GWAS..."
    df_gwas = pd.read_csv(file_gwas, index_col=False, header=0, delimiter="\s+") # read
    print "END: reading GWAS..."
    # rsID    beta    se      pval    snp_maf chr     pos
    # rs3107975       -0.0126339      0.00613983      0.095   0.01292 1       55326
    # rs114608975     0.00190558      0.0017629       0.27    0.06064 1       86028
    # rs28863004      -0.00409838     0.00717661      0.72    0.01789 1       526736
    # rs1988726       -0.0219751      0.0151771       0.2     0.04473 1       564862
    # rs371431021     -0.00330813     0.00776404      0.68    0.01889 1       565130

    n_snps_not_found_in_hapmap = sum(~df_gwas["rsID"].isin(df_collection["rsID"]))
    print "Number of SNPs in GWAS data: {}".format(len(df_gwas))
    print "Number of SNPs in GWAS data *NOT found* in HAPMAP file: {}".format(n_snps_not_found_in_hapmap)
    print "Percent SNPs not found: {:.2f} %".format(n_snps_not_found_in_hapmap/float(len(df_gwas))*100)


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
    file_out = "{}.gwassumstats.hapmap3.rolypoly_fmt.tab.gz".format(file_out_prefix)
    df_join.to_csv(file_out, sep="\t", index=False, compression='gzip')
    print "END: exported file: {}".format(file_out)

    list_res = [file_out_prefix,
                len(df_gwas), 
                n_snps_not_found_in_hapmap, 
                n_snps_not_found_in_hapmap/float(len(df_gwas))*100, 
                df_join.shape]
    return list_res


def main_parallel(file_gwas):
    file_out_prefix = os.path.basename(file_gwas).split(".")[0] # e.g. disease_RESPIRATORY_ENT.sumstats.gz --> disease_RESPIRATORY_ENT
    print "RUNNING {}".format(file_out_prefix)
    list_res = main(file_gwas, file_out_prefix)
    return list_res


def main_parallel2(i):
    file_gwas = list_file_gwas[i]
    file_out_prefix = list_file_output_prefix[i]
    print "RUNNING {}".format(file_out_prefix)
    list_res = main(file_gwas, file_out_prefix)
    return list_res



#####################################################################################
######################################## XXXX  ######################################
#####################################################################################

# list_gwas = ["body_HEIGHT_Wood2014",
#              "body_BMI_Yengo2018",
#              "body_height_Yengo2018"]

# list_files = ["/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/{}.gwassumstats.rolypoly_fmt.tab.gz".format(x) for x in list_gwas]


list_files = glob.glob("/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/*.gwassumstats.rolypoly_fmt.tab.gz")

print list_files

# list_file_output_prefix = ["NULL_GWAS_1KG_phase3_EUR_N{}".format(i) for i in range(1,10+1)]


#####################################################################################
######################################## MAIN  ######################################
#####################################################################################

df_collection = read_hapmap_snplist()

### Normal mode
# main(file_gwas, file_out_prefix)

### Alkes/loop mode
# for file_gwas in list_files:
#     file_out_prefix = os.path.basename(file_gwas).split(".")[0]
#     print "RUNNING {}".format(file_out_prefix)
#     main(file_gwas, file_out_prefix)

### Alkes/loop mode - parallel
### REMEMBER TO MAKE list_files 
print "Starting pool..."
pool = multiprocessing.Pool(80)
list_res = pool.map(main_parallel, list_files)
print list_res
for list_x in list_res:
    print "\t".join(str(x) for x in list_x)

## Loop mode - parallel2 | USE WEHN list_file_output_prefix is DEFINED
# print "Starting pool..."
# pool = multiprocessing.Pool(10)
# pool.map(main_parallel2, range(len(list_file_gwas)))




print "SCRIPT ENDED"

