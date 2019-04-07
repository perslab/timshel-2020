#!/usr/bin/env bash

### DEPICT must run on python. By loading python2 for this bash script we ensure that all its subprocesses inherit the correct python version/environment:
module unload anaconda
module load anaconda/2-4.4.0

### CONSTANTS
PY_SCRIPT=run_depict_from_template_cfg.py
CFG_TEMPLATE=/projects/timshel/sc-genetics/sc-genetics/src/depict/TEMPLATE_LDSC_MUNGED_GWAS.cfg


###################################### mousebrain_sem_mean ######################################

EXPRESSION_FILE=/projects/timshel/sc-genetics/sc-genetics/data/genes_cell_type_specific/mousebrain_all.mean.tab
EXPRESSION_NAME=mousebrain_sem_mean

GWAS=BMI_UPDATE_Yengo2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

GWAS=BMI_UKBB_Loh2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

wait

###################################### tabula_muris_sem_mean ######################################

EXPRESSION_FILE=/projects/timshel/sc-genetics/sc-genetics/data/genes_cell_type_specific/tabula_muris.mean.tab
EXPRESSION_NAME=tabula_muris_sem_mean

GWAS=BMI_UPDATE_Yengo2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

GWAS=BMI_UKBB_Loh2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

wait


###################################### STANDARD DEPICT TISSUE ######################################

EXPRESSION_FILE=/scratch/tmp-depict-1.0.174/data/tissue_expression/GPL570EnsemblGeneExpressionPerTissue_DEPICT20130820_z_withmeshheader.txt
EXPRESSION_NAME=depict_tissues

GWAS=BMI_UPDATE_Yengo2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

GWAS=BMI_UKBB_Loh2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

wait

###################################### tabula_muris.tissue_celltype.depict.z_score_two_step ######################################

EXPRESSION_FILE=/projects/timshel/DEPICT/data-expression/tabula_muris.tissue_celltype.depict.z_score_two_step.tab
EXPRESSION_NAME=tabula_muris_z_score_two_step

GWAS=BMI_UPDATE_Yengo2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

GWAS=BMI_UKBB_Loh2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

wait

###################################### mousebrain.all.depict.z_score_two_step ######################################

EXPRESSION_FILE=/projects/timshel/DEPICT/data-expression/mousebrain.all.depict.z_score_two_step.tab
EXPRESSION_NAME=mousebrain_z_score_two_step

GWAS=BMI_UPDATE_Yengo2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

GWAS=BMI_UKBB_Loh2018_no_mhc
python ${PY_SCRIPT} --cfg_template_file $CFG_TEMPLATE \
--gwas_summary_statistics_file /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/$GWAS.sumstats.gz \
--label_for_output_files $GWAS.5e-8.$EXPRESSION_NAME \
--tissue_expression_file $EXPRESSION_FILE &

wait


