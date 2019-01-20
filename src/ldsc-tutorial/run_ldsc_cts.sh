#!/bin/bash

# EXAMPLE:
# time run_ldsc_cts.sh PASS_BMI1.sumstats GTEx BMI1.GTEx
# time run_ldsc_cts.sh PASS_BMI1.sumstats GTEx BMI1.GTEx_no_baseline

sumstats=$1 #LDSC formatted sumstats file
cts_name=$2 #options: Cahoy, Franke, GTEx, GTEx_brain, [ImmGen, Roadmap,]
out=$3 # e.g. 

PATH_LDSC_DATA="/projects/timshel/sc-genetics/ldsc/data" # no trailing slash
PATH_LDSC_SEG="/projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores" # no trailing slash
PATH_GWAS="/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/alkesgroup/sumstats_formatted" # no trailing slash
PATH_LDSC_SRC="/projects/timshel/sc-genetics/ldsc/ldsc"

# run LDSC
$PATH_LDSC_SRC/ldsc.py \
    --h2-cts $PATH_GWAS/$sumstats \
    --ref-ld-chr $PATH_LDSC_SEG/$cts_name/$cts_name.control.,$PATH_LDSC_DATA/1000G_EUR_Phase3_baseline/baseline. \
    --out $out \
    --ref-ld-chr-cts $PATH_LDSC_SEG/$cts_name.ldcts \
    --w-ld-chr $PATH_LDSC_DATA/weights_hm3_no_hla/weights.

# # run LDSC - NO BASELINE
# $PATH_LDSC_SRC/ldsc.py \
#     --h2-cts $PATH_GWAS/$sumstats \
#     --ref-ld-chr $PATH_LDSC_SEG/$cts_name/$cts_name.control. \
#     --out $out \
#     --ref-ld-chr-cts $PATH_LDSC_SEG/$cts_name.ldcts \
#     --w-ld-chr $PATH_LDSC_DATA/weights_hm3_no_hla/weights.