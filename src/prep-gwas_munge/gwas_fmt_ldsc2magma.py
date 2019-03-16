#!/usr/bin/env python2.7

import os
import sys
import numpy as np
import pandas as pd
import time

import glob

import multiprocessing

import pdb

from scipy import stats # for stats.norm.cdf

###################################### USAGE ######################################

# time python gwas_to_rolypoly_fmt-map_chrpos.py |& tee log.alkes_group_UKBB.gwas_to_rolypoly_fmt-map_chrpos.txt


###################################### SCRIPT ######################################

def main(file_gwas, file_out_prefix):

    ### Read GWAS data
    print "START: reading GWAS..."
    df_gwas = pd.read_csv(file_gwas, index_col=False, header=0, delimiter="\s+") # read
    print "END: reading GWAS..."
    # SNP     A1      A2      N       CHISQ   Z
    # rs6671356       C       T       70100.0 0.172612905312  0.415467092935
    # rs6604968       G       A       70100.0 0.291125788806  0.539560736902
    # rs4970405       A       G       70100.0 0.102204513891  0.319694407037
    # rs12726255      G       A       70100.0 0.312418295691  0.558943911042
    # rs4970409       G       A       70100.0 0.0524226849517 0.228960007319


    ### Drop NA
    n_snps_start = df_gwas.shape[0]
    df_gwas.dropna(axis=0, how='any', inplace=True) # do this before converting "N" column to integer to avoid error "ValueError: Cannot convert non-finite values (NA or inf) to integer"
    n_snps_after_remove_na_rows = df_gwas.shape[0]
    n_na_snps_removed = n_snps_start - n_snps_after_remove_na_rows
    
    ### Convert N to integer column
    df_gwas["N"] = df_gwas["N"].astype(int)
    
    ### Get duplicates rsIDs
    bool_duplicated_snps = df_gwas.duplicated(subset="SNP", keep=False) # returns bool
    # keep : {'first', 'last', False}, default 'first'
    # first : Mark duplicates as True except for the first occurrence.
    # last : Mark duplicates as True except for the last occurrence.
    # False : Mark all duplicates as True.
    n_duplicated_snps = bool_duplicated_snps.sum()
    
    ### Calc P-value
    df_gwas["P"] = stats.norm.cdf(df_gwas["Z"]) # Python calculates left/lower-tail probabilities by default.
    # *OBS*: perhaps you should get the DOUBLE TAILED P-value instead! 
    
    ### Drop dups and select cols
    df_gwas_clean = df_gwas.loc[~bool_duplicated_snps, ["SNP", "P", "N"]]
    
    ### Drop NA again - just to be sure
    df_gwas_clean.dropna(axis=0, how='any', inplace=True) # do this before converting "N" column to integer to avoid error "ValueError: Cannot convert non-finite values (NA or inf) to integer"
    n_snps_final = df_gwas_clean.shape[0]


    print "Number of SNPs in GWAS data START: {}".format(n_snps_start)
    print "Number of SNPs in GWAS data FINAL: {}".format(n_snps_final)
    print "Number of duplicated SNPs removed: {}".format(n_duplicated_snps)
    print "Number of NA SNPs removed: {}".format(n_na_snps_removed)


    ## Export
    print "START: exporting file..."
    file_out = "{}.magma_fmt.txt.gz".format(file_out_prefix)
    df_gwas_clean.to_csv(file_out, sep="\t", index=False, compression='gzip')
    print "END: exported file: {}".format(file_out)

    list_res = [file_out_prefix,
                n_snps_start,
                n_snps_final,
                n_duplicated_snps,
                n_na_snps_removed]
    return list_res


def main_parallel(file_gwas):
    file_out_prefix = os.path.basename(file_gwas).split(".")[0] # e.g. disease_RESPIRATORY_ENT.sumstats.gz --> disease_RESPIRATORY_ENT
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

list_files = glob.glob("/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/finucane-collection/public_sumstats/*.sumstats")

print list_files



#####################################################################################
######################################## MAIN  ######################################
#####################################################################################


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





print "SCRIPT ENDED"

