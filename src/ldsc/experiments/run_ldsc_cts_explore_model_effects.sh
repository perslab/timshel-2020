

###################################### BMI on GTEx: effect of all_genes and baseline ######################################


### nothing, but here you need to loop over each annotation # SEE wrapper_run_ldsc_h2.py
# time /tools/anaconda/2-4.4.0/bin/python2 /projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py \
#     --h2 /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018.sumstats.gz \
#     --ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.1. \
#     --out /projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment.gtex_tissues.tissueX.no_baseline.no_all_genes.BMI_Yengo2018 \
#     --w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. \
#     --frqfile-chr /projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_frq/1000G.EUR.QC. \
#     --overlap-annot \
#     --print-coefficients



### baseline (no all genes)
time /tools/anaconda/2-4.4.0/bin/python2 /projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py \
    --h2-cts /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018.sumstats.gz \
    --ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline. \
    --out /projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment.gtex_tissues.baseline.BMI_Yengo2018 \
    --ref-ld-chr-cts /projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/ldsc.gtex_tissues_without_all_genes.ldscts.txt \
    --w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. &


### all genes (no baseline)
time /tools/anaconda/2-4.4.0/bin/python2 /projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py \
    --h2-cts /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018.sumstats.gz \
    --ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control. \
    --out /projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment.gtex_tissues.all_genes.BMI_Yengo2018 \
    --ref-ld-chr-cts /projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/ldsc.gtex_tissues_with_all_genes.ldscts.txt \
    --w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. &
# ---> produceses Singular matrix in regression. Don't know why.
# Reading cts reference panel LD Score from /projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.1.,/projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control.[1-22] ...
# Performing regression.
# Traceback (most recent call last):
#   File "/projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py", line 646, in <module>
#     sumstats.cell_type_specific(args, log)
#   File "/projects/timshel/sc-genetics/ldsc/ldsc/ldscore/sumstats.py", line 300, in cell_type_specific
#     twostep=None, old_weights=True)
#   File "/projects/timshel/sc-genetics/ldsc/ldsc/ldscore/regressions.py", line 346, in __init__
#     slow=slow, step1_ii=step1_ii, old_weights=old_weights)
#   File "/projects/timshel/sc-genetics/ldsc/ldsc/ldscore/regressions.py", line 208, in __init__
#     jknife = jk.LstsqJackknifeFast(x, y, n_blocks)
#   File "/projects/timshel/sc-genetics/ldsc/ldsc/ldscore/jackknife.py", line 309, in __init__
#     self.est = self.block_values_to_est(xty, xtx)
#   File "/projects/timshel/sc-genetics/ldsc/ldsc/ldscore/jackknife.py", line 386, in block_values_to_est
#     return np.linalg.solve(xtx, xty).reshape((1, p))
#   File "/tools/anaconda/2-4.4.0/lib/python2.7/site-packages/numpy/linalg/linalg.py", line 390, in solve
#     r = gufunc(a, b, signature=signature, extobj=extobj)
#   File "/tools/anaconda/2-4.4.0/lib/python2.7/site-packages/numpy/linalg/linalg.py", line 89, in _raise_linalgerror_singular
#     raise LinAlgError("Singular matrix")
# LinAlgError: Singular matrix


### baseline + all genes
time /tools/anaconda/2-4.4.0/bin/python2 /projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py \
    --h2-cts /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018.sumstats.gz \
    --ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline. \
    --out /projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment.gtex_tissues.baseline.all_genes.BMI_Yengo2018 \
    --ref-ld-chr-cts /projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/ldsc.gtex_tissues_with_all_genes.ldscts.txt \
    --w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. &



###################################### SCZ with and without liability ######################################
# ### CONCLUSION: no effect on taus of adding "--pop-prev 0.01 --samp-prev 0.43" !

# ### with liability (--pop-prev 0.01 --samp-prev 0.43)
# GWAS=SCZ_Ripke2014
# CTS=wgcna.nn_lira_sema
# time /tools/anaconda/2-4.4.0/bin/python2 /projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py \
#     --h2-cts /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}.sumstats.gz \
#     --ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline. \
#     --out /projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/${CTS}.SCZ_Ripke2014.liability \
#     --ref-ld-chr-cts /projects/timshel/sc-genetics/sc-genetics/src/ldsc/${CTS}.ldcts.txt \
#     --w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. \
#     --pop-prev 0.01 --samp-prev 0.43 &

# ### without liability
# GWAS=SCZ_Ripke2014
# CTS=wgcna.nn_lira_sema
# time /tools/anaconda/2-4.4.0/bin/python2 /projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py \
#     --h2-cts /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${GWAS}.sumstats.gz \
#     --ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline. \
#     --out /projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/${CTS}.SCZ_Ripke2014.no_liability \
#     --ref-ld-chr-cts /projects/timshel/sc-genetics/sc-genetics/src/ldsc/${CTS}.ldcts.txt \
#     --w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. &


# wait
