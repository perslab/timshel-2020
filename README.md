# Mapping heritability of obesity by cell types
Code to reproduce analysis and figures in [Timshel (2019, bioRxiv)](https://www.biorxiv.org/XXXX) **[TODO: update link]**: _Mapping heritability of obesity by cell types_.


### How to reproduce analysis
1. Download and prepare single-cell expression datasets: `src/prep-expr_datasets/<dataset>-preprocessing.R`
2. Perform CELLEX-R* ES precomputation: `src/es/run_es_precalculation.ipynb`
3. Calculate CELLEX-R* ES: `src/es/run_es_calculation.sh`
4. Run [CELLECT](https://github.com/perslab/CELLECT) analysis:

   CELLECT-LDSC config files in `src/CELLECT-LDSC`  
   CELLECT-MAGMA using scripts in `src/CELLECT-MAGMA`
5. Generate figures and tables using scripts in `src/publication`

* Since publication, we have implemented CELLEX-PY available [here](https://github.com/perslab/CELLEX). We recommend using CELLEX-PY as it is easier to use for new users.

Code to reproduce cell-type co-expression network analysis can be found here: https://github.com/perslab/19-BMI-brain-wgcna
