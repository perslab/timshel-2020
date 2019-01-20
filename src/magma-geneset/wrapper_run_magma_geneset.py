

import sys
import subprocess
import os
import glob


### Usage
# python run_magma_geneset.py \
# --file_gwas_raw /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw \
# --file_wgcna_module_genelist /raid5/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv \
# --outprefix magma_geneset_BMI_Yengo2018

# ll /scratch/tmp-magma_gwas/*.genes.raw

list_gwas = ["BMI_Yengo2018",
"WHR_adjBMI_Shungin2015",
"WHR_Shungin2015",
"EA3_Lee2018",
"HEIGHT_Wood2014",
"SCZ_Ripke2014",
"RA_Okada2014",
"1KG_phase3_EUR_null_gwas_P10"]

# list_gwas = ["RA_Okada2014"]



for gwas in list_gwas:
    print(gwas)
    ### Novo bulk
    # cmd = "python run_magma_geneset.py --file_gwas_raw /scratch/tmp-magma_gwas/{gwas}.txt.10UP.1.5DOWN.genes.raw --file_wgcna_module_genelist /raid5/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv --outprefix magma_geneset_{gwas}".format(gwas=gwas)
    ### Mousebrain neurons
    # cmd = "python run_magma_geneset.py --file_gwas_raw /scratch/tmp-magma_gwas/{gwas}.txt.10UP.1.5DOWN.genes.raw --file_wgcna_module_genelist /projects/jonatan/tmp-mousebrain/tables/mousebrain_Neurons_ClusterName_2_cell_cluster_module_genes.csv --outprefix magma_geneset_{gwas}".format(gwas=gwas)
    ### MACA
    cmd = "python run_magma_geneset.py --file_gwas_raw /scratch/tmp-magma_gwas/{gwas}.txt.10UP.1.5DOWN.genes.raw --file_wgcna_module_genelist /projects/jonatan/tmp-maca/tables/maca_tissue_cell_type_kME_cell_cluster_module_genes.csv --outprefix magma_geneset_{gwas}".format(gwas=gwas)
    subprocess.Popen(cmd, shell=True)

print("WRAPPER SCRIPT DONE")