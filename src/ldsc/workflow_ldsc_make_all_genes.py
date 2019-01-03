
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import re


##################################################################################################
###################################### PARAMS AND CONSTANTS ######################################
##################################################################################################


PYTHON_EXEC = "/tools/anaconda/2-4.4.0/bin/python2" # runs on python2
PATH_LDSC_SCRIPT = "/raid5/projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py" 


################## CONTROL GENE SET ##################
FLAG_BINARY = True # does not matter, but let's just keep the numbers binary
FLAG_WGCNA = False

### FDR significant modules
prefix_genomic_annot="control.all_genes_in_dataset"
file_multi_gene_set="/raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.all_genes_in_dataset.txt"
### SNIPPET [the annotation name does NOT matter when the ldscores are only used for control [but the name will be used when doing --h2 (non cts analysis)]]
# all_genes_in_dataset.mousebrain      ENSG00000141668 1
# all_genes_in_dataset.mousebrain      ENSG00000204624 1
# all_genes_in_dataset.mousebrain      ENSG00000187848 1
# all_genes_in_dataset.mousebrain      ENSG00000171522 1


#########################################################################################
###################################### PRE-PROCESS ######################################
#########################################################################################

### Make annot
###  *RESOURCE NOTES*: if you have many modules (~3000-6000) then set --n_parallel_jobs to ~2-5 (instead of 22). Otherwise the script will up all the MEMORY on yggdrasil and fail.
cmd = """{PYTHON_EXEC} make_annot_from_geneset_all_chr.py \
--file_multi_gene_set {file_multi_gene_set} \
--file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
--windowsize 100000 \
--bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
{flag_binary} \
{flag_wgcna} \
--out_dir /scratch/sc-ldsc/{prefix_genomic_annot} \
--out_prefix {prefix_genomic_annot}
""".format(PYTHON_EXEC=PYTHON_EXEC, 
	file_multi_gene_set=file_multi_gene_set, 
	prefix_genomic_annot=prefix_genomic_annot, 
	flag_wgcna="--flag_wgcna --flag_mouse" if FLAG_WGCNA else "",
	flag_binary="--flag_encode_as_binary_annotation" if FLAG_BINARY else "",
	) 
# --n_parallel_jobs 11

print("Running command: {}".format(cmd))
p = subprocess.Popen(cmd, shell=True)
p.wait()
print("Return code: {}".format(p.returncode))
if not p.returncode == 0:
	raise Exception("Got non zero return code")


### compute LD scores
### *RESOURCE NOTES*: this script uses a lot of CPU. Never run more than 4 parallel jobs. 4 parallel jobs will use ~220% CPU
cmd="{PYTHON_EXEC} wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/{prefix_genomic_annot}/ --n_parallel_jobs 2".format(PYTHON_EXEC=PYTHON_EXEC, prefix_genomic_annot=prefix_genomic_annot)
print("Running command: {}".format(cmd))
p = subprocess.Popen(cmd, shell=True)
p.wait()
print("Return code: {}".format(p.returncode))
# RUNTIME ----> ~6 h for ~500 modules with --n_parallel_jobs=4
if not p.returncode == 0:
	raise Exception("Got non zero return code")

### split LD scores
### This script will read 1 ".COMBINED_ANNOT.$CHR.l2.ldscore.gz" file  (N_SNPs x N_Modules) per parallel process.
###  *RESOURCE NOTES*: this script may use quiet a lot of memory for many modules? Not sure
cmd="{PYTHON_EXEC} split_ldscores.py --prefix_ldscore_files /scratch/sc-ldsc/{prefix_genomic_annot}/ --n_parallel_jobs 4".format(PYTHON_EXEC=PYTHON_EXEC, prefix_genomic_annot=prefix_genomic_annot)
print("Running command: {}".format(cmd))
p = subprocess.Popen(cmd, shell=True)
p.wait()
print("Return code: {}".format(p.returncode))
# RUNTIME ----> ~10 min
if not p.returncode == 0:
	raise Exception("Got non zero return code")

###################################### XXXX ######################################

print("Script is done!")



