# Mapping heritability of obesity by cell types

This repository contains code to reproduce analysis and figures in presented in:  

**[Timshel (bioRxiv, 2020): _Mapping heritability of obesity by cell types_](https://www.biorxiv.org/content/10.1101/2020.01.27.920033v1)**

![image](https://user-images.githubusercontent.com/5487016/72666147-aecfdd00-3a0f-11ea-8609-12c157ee24f3.png)

## Steps to reproduce analysis
1. Download and prepare single-cell expression datasets: `src/prep-expr_datasets/<dataset>-preprocessing.R`
2. Perform [CELLEX](https://github.com/perslab/CELLEX)* ES precomputation: `src/es/run_es_precalculation.ipynb`
3. Calculate CELLEX* ES: `src/es/run_es_calculation.sh`
4. Run [CELLECT](https://github.com/perslab/CELLECT) analysis:

   CELLECT-LDSC config files in `src/CELLECT-LDSC`  
   CELLECT-MAGMA using scripts in `src/CELLECT-MAGMA`  

   Note that to run these analysis, you will have to download the relevant GWAS summary stats from the consortium hosting the files. We are not allowed to share (redistribute) the the summary statistics. Details on the GWAS summary stats and how they where 'munged' can be found in `src/prep-gwas_munge/README-cmds_munge_sumstats.txt`
5. Generate figures and tables using the scripts in [`src/publication`](https://github.com/perslab/timshel-bmicelltypes/tree/master/src/publication)**


\*Since finishing the manuscript, we have implemented CELLEX in Python available [here](https://github.com/perslab/CELLEX). We recommend using Python CELLEX as it is easier to use for new users.

\*\*For cross-platform portability we used [renv](https://rstudio.github.io/renv/articles/renv.html) to create a reproducible environment for the R code used to produce the figures/tables. You can install the environment (i.e. CRAN packages) by running either `renv::init()` or `renv::restore(lockfile="./renv.lock")`. (See the renv.lock file [here](https://github.com/perslab/timshel-bmicelltypes/blob/master/renv.lock).)


## Additional repositories to reproduce analysis

Code to reproduce cell-type co-expression networks can be found here: [19-BMI-brain-wgcna](https://github.com/perslab/19-BMI-brain-wgcna)

Code to reproduce 'obesity geneset' enrichment tests can be found here: [19-BMI-brain-genesettests](https://github.com/perslab/19-BMI-brain-genesettests)

![image](https://user-images.githubusercontent.com/5487016/72666162-d9219a80-3a0f-11ea-94c2-669125fd588a.png)

## Figures and tables

All output figures (.pdf) and tables (.csv) are available in this repository: [src/publication/figs](https://github.com/perslab/timshel-bmicelltypes/tree/master/src/publication/figs) and [src/publication/tables](https://github.com/perslab/timshel-bmicelltypes/tree/master/src/publication/tables).





