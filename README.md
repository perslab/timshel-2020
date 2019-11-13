# Mapping heritability of obesity by cell types

This repository contains code to reproduce analysis and figures in presented in:  

**[Timshel (bioRxiv, 2019)](https://www.biorxiv.org/XXXX): _Mapping heritability of obesity by cell types_.**


### Steps to reproduce analysis
1. Download and prepare single-cell expression datasets: `src/prep-expr_datasets/<dataset>-preprocessing.R`
2. Perform CELLEX-R* ES precomputation: `src/es/run_es_precalculation.ipynb`
3. Calculate CELLEX-R* ES: `src/es/run_es_calculation.sh`
4. Run [CELLECT](https://github.com/perslab/CELLECT) analysis:

   CELLECT-LDSC config files in `src/CELLECT-LDSC`  
   CELLECT-MAGMA using scripts in `src/CELLECT-MAGMA`  

   Note that to run these analysis, you will have to download the relevant GWAS summary stats from the consortium hosting the files. We are not allowed to share (redistribute) the the summary statistics. Details on the GWAS summary stats and how they where 'munged' can be found in `src/prep-gwas_munge/README-cmds_munge_sumstats.txt`
5. Generate figures and tables using scripts in `src/publication`


\*Since publication, we have implemented CELLEX-PY available [here](https://github.com/perslab/CELLEX). We recommend using CELLEX-PY as it is easier to use for new users.

### Additional repositories to reproduce analysis

Code to reproduce cell-type co-expression networks can be found here: https://github.com/perslab/19-BMI-brain-wgcna

Code to reproduce 'obesity geneset' enrichment tests can be found here: https://github.com/perslab/19-BMI-brain-genesettests
