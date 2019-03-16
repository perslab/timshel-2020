
import os
import sys
import numpy as np
import pandas as pd
import time

import glob


################## alkes_lo_UKBB ##################

# SNP            CHR      POS     A1      A2      REF     EAF             Beta            se              P       N       INFO
# rs10399793      1       49298   T       C       T       0.37622         0.000334086     0.000731539     7.0E-01 459324  0.342797
# rs2462492       1       54676   C       T       C       0.599409        -0.000937692    0.000724671     1.7E-01 459324  0.340158
# rs3107975       1       55326   T       C       T       0.991552        0.00706826      0.0040343       8.3E-02 459324  0.324228
# rs74447903      1       57033   T       C       T       0.998221        0.00919789      0.00897642      3.8E-01 459324  0.296256
# 1:70728_C_T     1       70728   C       T       C       0.997834        0.00859618      0.00730365      2.6E-01 459324  0.365713
# rs2462495       1       79033   A       G       A       0.00129115      0.00321513      0.00929391      7.7E-01 459324  0.536566
# rs114608975     1       86028   T       C       T       0.896384        0.000549046     0.00115835      5.8E-01 459324  0.340885

list_files = glob.glob("/Users/djw472/data/GWAS-sumstats/alkesgroup-collection/UKBB/*sumstats.gz")
print list_files
    # /Users/djw472/data/GWAS-sumstats/alkesgroup-collection/UKBB/body_BMIz.sumstats.gz
    # /Users/djw472/data/GWAS-sumstats/alkesgroup-collection/UKBB/body_HEIGHTz.sumstats.gz

COL_rsID = "SNP"
COL_beta = "Beta"
COL_se = "se"
COL_pval = "P" # this is not strictly needed

COL_chr = "CHR"
COL_pos = "POS"
COL_maf = "EAF"


###################################### FUNCTIONS ######################################


def main(file_collection, file_out_prefix):

    ### Read GWAS data
    print "START: reading GWAS..."
    df_gwas = pd.read_csv(file_gwas, usecols=GWAS_COL_SPECS.keys(), index_col=False, header=0, delimiter="\t") # read
    df_gwas = df_gwas.rename(columns=GWAS_COL_SPECS) # rename cols. REF: https://pandas.pydata.org/pandas-docs/stable/basics.html#basics-rename
    print "END: reading GWAS..."
    # df_gwas.head()
    #     rsID    beta    se  pval
    # 0   rs1000000   0.0001  0.0043  0.98140
    # 1   rs10000010  -0.0022 0.0029  0.43840
    # 2   rs10000012  -0.0096 0.0053  0.07009
    # 3   rs10000013  -0.0096 0.0043  0.02558
    # 4   rs10000017  -0.0038 0.0045  0.39840


    print "Number of SNPs in GWAS data: {}".format(len(df_gwas))

    ## Export
    print "START: exporting file..."
    file_out = "{}.gwassumstats.rolypoly_fmt.tab.gz".format(file_out_prefix)
    df_join.to_csv(file_out, sep="\t", index=False, compression='gzip')
    print "END: exported file: {}".format(file_out)




########################################################################################
###################################### GWAS PARAMETERS ######################################
########################################################################################

##################  ##################
### SCZ
# hg19chrc        snpid   a1      a2      bp      info    or      se      p       ngt
# chr1    rs4951859       C       G       729679  0.631   0.97853 0.0173  0.2083  0
# chr1    rs142557973     T       C       731718  0.665   1.01949 0.0198  0.3298  0
# chr1    rs141242758     T       C       734349  0.666   1.02071 0.02    0.3055  0


#####################################################################################
###################################### CONSTANTS ######################################
#####################################################################################


GWAS_COL_SPECS = {
    COL_rsID:"rsID",
    COL_beta:"beta",
    COL_se:"se",
    COL_pval:"pval",
    COL_chr:"chr",
    COL_pos:"pos",
    COL_maf:"maf"
}


#####################################################################################
######################################## MAIN  ######################################
#####################################################################################

## Alkes mode
for file_gwas in list_files:
    file_out_prefix = os.path.basename(file_gwas).split(".")[0]
    print "RUNNING {}".format(file_out_prefix)
    main(file_gwas, file_out_prefix)





print "SCRIPT ENDED"

