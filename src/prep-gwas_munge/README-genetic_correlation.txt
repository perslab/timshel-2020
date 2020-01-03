# Commands for genetic correlation between BMI and Hippoc

FILE1=/projects/timshel/data/gwas_sumstats_ldsc/timshel-collection/HIPPOVOL_Hilbar2017.sumstats.gz
FILE2=/projects/timshel/data/gwas_sumstats_ldsc/timshel-collection/BMI_UKBB_Loh2018.sumstats.gz
LDSC_DATA_DIR=/nfsdata/projects/timshel/sc-genetics/ldsc/data
LDSC_SRC_DIR=/nfsdata/projects/timshel/sc-genetics/CELLECT/ldsc
FILEOUT=/nfsdata/projects/timshel/sc-genetics/timshel-bmicelltypes2019/src/prep-gwas_munge/ldsc_out.BMI_UKBB_Loh2018__HIPPOVOL_Hilbar2017

python2 $LDSC_SRC_DIR/ldsc.py \
--rg $FILE1,$FILE2 \
--ref-ld-chr $LDSC_DATA_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSC_DATA_DIR/eur_w_ld_chr/ \
--out $FILEOUT


### Info from LDSC WIKI: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
 For partitioned LD Score regression, you will want to have separate sets of LD Scores for the
 --w-ld-chr and --ref-ld-chr flags. We also used different LD Scores for non-partitioned LD Score
 regression in the --w-ld-chr and --ref-ld-chr flags in the original LD Score regression paper ("LD
 Score regression distinguishes confounding from polygenicity in genome-wide association studies.")
 However, ***our current recommendation is to use the same LD Scores for both flags for non-partitioned
 LD Score regression. This is the procedure described in this tutorial***.


###################################### DID NOT WORK ######################################

---> The below did not work because 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. does not contain the "M_5_50" files.

FILE1=/projects/timshel/data/gwas_sumstats_ldsc/timshel-collection/HIPPOVOL_Hilbar2017.sumstats.gz
FILE2=/projects/timshel/data/gwas_sumstats_ldsc/timshel-collection/BMI_UKBB_Loh2018.sumstats.gz
LDSC_DATA_DIR=/nfsdata/projects/timshel/sc-genetics/CELLECT/data/ldsc
LDSC_SRC_DIR=/nfsdata/projects/timshel/sc-genetics/CELLECT/ldsc
FILEOUT=/nfsdata/projects/timshel/sc-genetics/timshel-bmicelltypes2019/src/prep-gwas_munge/ldsc_out.BMI_UKBB_Loh2018__HIPPOVOL_Hilbar2017

python2 $LDSC_SRC_DIR/ldsc.py \
--rg $FILE1,$FILE2 \
--ref-ld-chr $LDSC_DATA_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--w-ld-chr $LDSC_DATA_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--out $FILEOUT
