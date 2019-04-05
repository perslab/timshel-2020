
from __future__ import print_function
import argparse
import glob
import os
import subprocess

###################################### USAGE ######################################
# Compatibility: Python 2 and 3

# time python wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/nn_lira_sema-combined/nn_lira_sema --n_parallel_jobs 4
# time python wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/nn_lira_sema-combined-fixed_print_snps/ --n_parallel_jobs 4
# time python wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/maca/ --n_parallel_jobs 1
# time python wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/mousebrain/ --n_parallel_jobs 1
# time python wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/hypothalamus_mette_thesis/ --n_parallel_jobs 4

###################################### DESCRIPTION ######################################

### Purpose
# Calculate LD scores for annotation files matching the input filepath prefix

### Steps
# 1) get list annot files
# 2) check that LDscores for annot files have not already been calculated
# 2) run subprocess: calculate LDscore for each annot file (annotation + 1..22 chr)

### Output
# This script will call ldsc.py with the --l2 command. 
# The following output files will be written to the --prefix_annot_files:
# <OUT>.l2.ldscore.gz
# <OUT>.l2.M
# <OUT>.l2.M_5_50
# <OUT>.log

### Example --l2 call to ldsc.py to compute ld scores:
# REF --thin-annot: https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial#step-2-computing-ld-scores-with-an-annot-file
# python2 ldsc.py \
#     --l2 \
#     --bfile /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.22 \
#     --ld-wind-cm 1 \
#     --annot GTEx_Cortex.annot.gz \
#     --thin-annot \
#     --out GTEx_Cortex \
#     --print-snps /projects/timshel/sc-genetics/ldsc/data/hapmap3_snps/hm.22.snp

###################################### NOTES ######################################

# Potential optimization: minimize I/O by calculating LDscores for the same chromosome across multiple annotations (multi-annotation file, ala baseline model). 
# Requires merging of annotation files across annotations. This would not be compatible with the "-ref-ld-chr-cts" workflow.

#### How much time/memory does it take to estimate LD Scores? 
# REF: https://github.com/bulik/ldsc/wiki/FAQ
# A. Using a 1cM window, it takes about 1 hour and 1GB RAM on my 1.7GHz Macbook Air to compute LD Scores for chromosome 1 using the N=378 1000 Genomes Europeans. 
# The time and memory complexity are both linear in the number of samples and the number of SNPs. 
# ***Computing partitioned LD Scores with a large number of partitions increases the memory usage (~8GB for 50 categories) but has only minimal impact on the runtime.***


### CPU usage:  
# * "It looks like ldsc will use all cores available when computing LD scores. Is there some ways to limit the number of threads to use?"
# * "how to control the num of threads for cpu" https://groups.google.com/forum/#!msg/ldsc_users/Bnj7FFl5jlw/nrQ7yH-8BgAJ


###################################### CONSTANTS ######################################


python_exec = "/tools/anaconda/2-4.4.0/bin/python2" # runs on python2
path_ldsc_script = "/projects/timshel/sc-genetics/ldsc/ldsc-timshel/ldsc.py" # important to use this version of the script, that will not compute the 'Annotation Correlation Matrix' to the log file
prefix_file_plink_per_chr = "/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC" # e.g. full file name will be <1000G.EUR.QC>.22.<bed,bim,fam>
# dir_file_print_snps = "/projects/timshel/sc-genetics/ldsc/data/hapmap3_snps" # full name will be <dir_file_print_snps>/hm.22.snp the LD scores are calculated/reported only for HapMap3 SNPs (--print-snps flag)
file_print_snps = "/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/print_snps.txt" # THIS is the file used for the "--print-snps" for ALL the baseline models. The files hm.1.snp, hm.2.snp, ..., does NOT contain all SNPs.

###################################### MAIN ######################################

parser = argparse.ArgumentParser()
parser.add_argument('--prefix_annot_files', type=str, help="""File path prefix to LDSC annot files to calculate ld scores. 
	Will calculate ld scores for annot files matching <prefix_annot_files>*.annot.gz.
	Notice that if prefix_annot_files is a directory (with a trailing '/'), then all annot file in the directory will be calculated.
	LD score output files will be written to the directory of prefix_annot_files [os.path.dirname(prefix_annot_files)]. 
	Example with prefix as input: /scratch/sc-ldsc/ldsc_analysis_dir/annot_files_prefix
	Example with dir as input: /scratch/sc-ldsc/ldsc_analysis_dir/ [OBS the trailing slash is important!]""")
parser.add_argument('--n_parallel_jobs', type=int, default=1, help='Number of parallel jobs to run.')


args = parser.parse_args()


files_annot = glob.glob("{}*.annot.gz".format(args.prefix_annot_files)) # e.g. <SOMEPATH>/nn_lira_sema.yellow.22.annot.gz | output files from make_annot_from_geneset_all_chr.py follow the pattern <PREFIX>.<ANNOTATIONNAME>.<CHR>.annot.gz
print("Found {} annot files matching prefix_annot_files".format(len(files_annot)))

###################################### CHECK FOR EXISTING LD SCORES ######################################

files_to_run = []
for file_annot in files_annot:
	file_ldscore = "{}.l2.ldscore.gz".format(file_annot.strip(".annot.gz"))
	if not os.path.exists(file_ldscore):
		files_to_run.append(file_annot)
	else:
		print("Annotation file {} has l2.ldscore.gz file. Will not compute ld scores again".format(file_annot))
print("Will compute for {} files".format(len(files_to_run)))

###################################### CALL ######################################

list_cmds = []
for file_annot in files_to_run:
	chromosome = file_annot.strip(".annot.gz").split(".")[-1] # file_annot = <SOMEPATH>/nn_lira_sema.yellow.22.annot.gz --> chromosome = 22
	file_plink_per_chr = "{}.{}".format(prefix_file_plink_per_chr, chromosome)

	# We also now allow for "thin annot" files, which omit the CHR, BP, SNP and CM columns and only have data on the annotation itself. 
	# The software assumes but does not check that thin annot files have the same SNPs in the same order as the plink files you used to compute LD scores. 
	# These require you to use the --thin-annot flag, as described below.
	cmd = "{python_exec} {script} --l2 --bfile {bfile} --ld-wind-cm 1 --annot {file_annot} --thin-annot --out {out} --print-snps {file_print_snps}".format(
		python_exec=python_exec,
		script=path_ldsc_script,
		bfile=file_plink_per_chr, # Prefix to plink data. E.g. <PATH>/1000G.EUR.QC.22 for chromosome 22
		file_annot=file_annot, # .annot.gz file
		out=file_annot.strip(".annot.gz"), # OBS: remember that the '--out' argument to ldsc.py is the FILEPATH PREFIX.
		file_print_snps=file_print_snps
		# file_print_snps="{}/hm.{}.snp".format(dir_file_print_snps, chromosome) # e.g. <dir_file_print_snps>/hm.22.snp
		)
	list_cmds.append(cmd)


# You need to keep devnull open for the entire life of the Popen object, not just its construction. 
FNULL = open(os.devnull, 'w') # devnull filehandle does not need to be closed?

list_of_processes = []
batch = 1
for i, cmd in enumerate(list_cmds, start=1):
	print("batch = {} | i = {} | Running command: {}".format(batch, i, cmd))
	p = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	list_of_processes.append(p)
	print("batch = {} | i = {} | PIDs of running jobs (list_of_processes):".format(batch, i))
	print(" ".join([str(p.pid) for p in list_of_processes])) # print PIDs
	if i % args.n_parallel_jobs == 0: # JOB BATCH SIZE
		batch += 1
		for p in list_of_processes:
			print("=========== Waiting for process: {} ===========".format(p.pid))
			p.wait()
			print("Returncode = {}".format(p.returncode))
		list_of_processes = [] # 'reset' list

### wait for the rest for the rest of the processes
for p in list_of_processes:
	print("=========== Waiting for process: {} ===========".format(p.pid))
	p.wait()



print("Script is done!")



