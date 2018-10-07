#!/usr/bin/env python2.7
from __future__ import print_function

### PT adds
import glob
import os
import subprocess

###################################### USAGE ######################################
# Must run on python2.7??




###################################### Algorithm description ######################################

# 1) get list annot files
# 2) check that LDscores for annot files have not already been calculated
# 2) run subprocess: calculate LDscore for each annot file (annotation + 1..22 chr)

### Example --l2 call to ldsc.py to compute ld scores:
# REF --thin-annot: https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial#step-2-computing-ld-scores-with-an-annot-file
# python2 ldsc.py \
#     --l2 \
#     --bfile /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.22 \
#     --ld-wind-cm 1 \
#     --annot GTEx_Cortex.annot.gz \
#     --thin-annot \
#     --out GTEx_Cortex \
#     --print-snps /raid5/projects/timshel/sc-genetics/ldsc/data/hapmap3_snps/hm.22.snp

###################################### NOTES ######################################

# Potential optimization: minimize I/O by calculating LDscores for the same chromosome across multiple annotations (multi-annotation file, ala baseline model). 
# Requires merging of annotation files across annotations. This would not be compatible with the "-ref-ld-chr-cts" workflow.

### CPU usage:  
# * "It looks like ldsc will use all cores available when computing LD scores. Is there some ways to limit the number of threads to use?"
# * "how to control the num of threads for cpu" https://groups.google.com/forum/#!msg/ldsc_users/Bnj7FFl5jlw/nrQ7yH-8BgAJ

###################################### FILE SNIPPETS ######################################





###################################### MAIN ######################################

parser = argparse.ArgumentParser()
parser.add_argument('--dir_ldsc_files', type=str, help='Path to directory with LDSC annot files to calculate LDSC. This will also be used as the output path')
# parser.add_argument('--n_proc', type=int, help='Number of processes. Default 22 (the number of chromosomes)')


args = parser.parse_args()

dir_ldsc_files = args.dir_ldsc_files

files_annot = glob.glob("*.annot.gz") # e.g. 


python_exec = "/tools/anaconda/2-4.4.0/bin/python2" # runs on python2
path_ldsc_script = "/raid5/projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py"
file_bim_prefix = "/raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.22"
file_print_snps = "/raid5/projects/timshel/sc-genetics/ldsc/data/hapmap3_snps/hm.22.snp" # the LD scores are calculated/reported only for HapMap3 SNPs (--print-snps flag)


cmd = "{exec} {script} --l2 --bfile {bfile} --ld-wind-cm 1 --annot {file_annot} --thin-annot --out {out} --print-snps {file_print_snps}".format(
	exec=python_exec,
	script=path_ldsc_script,
	bfile=file_bim, ### <-- UPDATE TO XXX
	file_annot=X, ### <-- UPDATE TO XXX
	out=, ### <-- UPDATE TO XXX
	file_print_snps=file_print_snps
	)

print("Making call: {}".format(cmd))
# p = subprocess.Popen(cmd, shell=True)
# p.wait()
# print("LDSC call done. Returncode = {}".format(p.returncode))


ldsc.py \
    --l2 \
    --bfile /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.22 \
    --ld-wind-cm 1 \
    --annot GTEx_Cortex.annot.gz \
    --thin-annot \
    --out GTEx_Cortex \
    --print-snps /raid5/projects/timshel/sc-genetics/ldsc/data/hapmap3_snps/hm.22.snp


print("Script is done!")



