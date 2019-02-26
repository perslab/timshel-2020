
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import re
import sys


import collections


###################################### USAGE ######################################
# Compatibility: Python 2 and 3

### Run in unbuffered mode
# time python -u workflow_ldsc_cts.py |& tee workflow_ldsc_cts.UNNAMED.out.txt

###################################### DOCS ######################################

### AIM: run LDSC --h2 with baseline model to estimate h2 observed. We will later use h2 estimate to calculate tau* and normalized-tau.


###################################### Job scheduler ######################################

def job_scheduler(list_cmds, n_parallel_jobs):
	""" Schedule parallel jobs with at most n_parallel_jobs parallel jobs."""
	list_of_processes = []
	batch = 1
	for i, cmd in enumerate(list_cmds, start=1):
		print("job schedule batch = {} | i = {} | Running command: {}".format(batch, i, cmd))
		## p = subprocess.Popen(cmd, shell=True, bufsize=0 if FLAG_UNBUFFERED else -1, stdout=FNULL, stderr=subprocess.STDOUT)
		### You need to keep devnull open for the entire life of the Popen object, not just its construction. 
		### FNULL = open(os.devnull, 'w') # devnull filehandle does not need to be closed?
		p = subprocess.Popen(cmd, shell=True, bufsize=0 if FLAG_UNBUFFERED else -1)
		list_of_processes.append(p)
		print("job schedule batch = {} | i = {} | PIDs of running jobs (list_of_processes):".format(batch, i))
		print(" ".join([str(p.pid) for p in list_of_processes])) # print PIDs
		if i % n_parallel_jobs == 0: # JOB BATCH SIZE
			batch += 1
			for p in list_of_processes:
				print("=========== Waiting for process: {} ===========".format(p.pid))
				sys.stdout.flush()
				p.wait()
				print("Returncode = {}".format(p.returncode))
			list_of_processes = [] # 'reset' list

	### wait for the rest for the rest of the processes
	for p in list_of_processes:
		print("=========== Waiting for process: {} ===========".format(p.pid))
		p.wait()

	return list_of_processes

##################################################################################################
############################################ CONSTANTS ###########################################
##################################################################################################



PYTHON3_EXEC = "/tools/anaconda/3-4.4.0/envs/py3_anaconda3_PT180510/bin/python3"
PYTHON2_EXEC = "/tools/anaconda/3-4.4.0/envs/py27_anaconda3_PT170705/bin/python2"

PATH_LDSC_SCRIPT = "/raid5/projects/timshel/sc-genetics/ldsc/ldsc-timshel/ldsc.py" 
PATH_LDSC_DATA_MAIN="/raid5/projects/timshel/sc-genetics/ldsc/data"


FLAG_UNBUFFERED = True

N_PARALLEL_LDSC_REGRESSION_JOBS = 50



##################################################################################################
############################################ PARAMS ##############################################
##################################################################################################


### Get all GWAS
list_gwas_filepaths = glob.glob("/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/*.sumstats.gz")
# dict_gwas = {re.search(r"(.*?)\.sumstats.gz", os.path.basename(filepath)).group(1):filepath for filepath in list_gwas_filepaths} # dict comprehension. e.g. BMI_UKBB_Loh2018:/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_UKBB_Loh2018.sumstats.gz

dict_gwas = {}
for filepath in list_gwas_filepaths:
	gwas = re.search(r"(.*?)\.sumstats.gz", os.path.basename(filepath)).group(1) # e.g. BMI_UKBB_Loh2018
	dict_gwas[gwas] = {"filepath": filepath}


######################################  ######################################
dict_gwas["SCZ_Pardinas2018_liability_scale"] = {"filepath": "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/SCZ_Pardinas2018.sumstats.gz",
								 "pop-prev":0.01, # from Gandal
								 "samp-prev":0.39} # CASES = 40,675 | CONTROLS = 64,643 | 40675/105318=.386211284

######################################  ######################################
dict_gwas["MS_Patsopoulos2011_liability_scale"] = {"filepath": "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/MS_Patsopoulos2011.sumstats.gz",
								 "pop-prev":0.001, # North America and Europe (>100/100,000 inhabitants) REF https://www.ncbi.nlm.nih.gov/pubmed/26718593
								 "samp-prev":0.31} # 5545/17698 | 5545 cases and 12153 controls | REF https://www.ncbi.nlm.nih.gov/pubmed/22190364 (MS_Patsopoulos2011)

######################################  ######################################
dict_gwas["RA_Okada201_liability_scale"] = {"filepath": "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/RA_Okada2014.sumstats.gz",
								 "pop-prev":0.01, # "The occurrence of RA is relatively constant with a prevalence of between 0.5 and 1.0%, a frequency that has been reported from several European [1-8] and North-American populations [9,10]" REF https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3240153/
								 "samp-prev":0.29} # 14361/58284=.246396953 | Eurpean RA GWAS meta-analysis (14361 RA cases and 43923 conrols) | N_total = 58284 | REF http://plaza.umin.ac.jp/~yokada/datasource/software.htm


######################################  ######################################
dict_gwas["INSOMNIA_Jansen2018_liability_scale"] = {"filepath": "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/INSOMNIA_Jansen2018.sumstats.gz",
								 "pop-prev":0.20, # hard to estimate exactly
								 "samp-prev":0.28} # The prevalence of insomnia was 28.3% in the UKB version 2 sample --> the gwas summary stats used is for UKBB only [N=386533].

### pop-prev
# REF: https://www.nature.com/articles/s41588-018-0333-3 (Jansen2018) | The diagnostic criteria for insomnia disorder2 (that is, difficulties with initiating or maintaining sleep with accompanying daytime complaints at least three times a week for at least three months, which cannot be attributed to inadequate circumstances for sleep3) are met by 10% of individuals, and up to one-third of older age individuals4.
# REF: https://www.nature.com/articles/s41380-018-0033-5 | Insomnia is highly prevalent, affecting 10–20% of adults in the United States [1] and worldwide [2].

### samp-prev | REF Jansen2018
# The UKB assessed insomnia complaints (hereafter referred to as ‘insomnia’)
# with a touchscreen device, whereas 23andMe research participants completed
# online surveys (Supplementary Tables 1 and 2). The assessment of insomnia in
# both samples shows high accuracy for insomnia disorder in the UKB and somewhat
# lower accuracy in 23andMe (sensitivity/specificity: UKB = 98/96%; 23andMe =
# 84/80%) (see Supplementary Note). The prevalence of insomnia was 28.3% in the
# UKB version 2 sample, 30.5% in the 23andMe sample, and 29.9% in the combined
# sample, which is in keeping with previous estimates for people of advanced age
# in the UK4 and elsewhere13,14.
######################################  ######################################




#########################################################################################
###################################### PRE-PROCESS ######################################
#########################################################################################

# ....

#########################################################################################
###################################### RUN LDSC h2 ######################################
#########################################################################################


### Create job commands
list_cmds_ldsc_prim = []
for gwas in dict_gwas:
	fileout_prefix_ldsc_h2 = "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2_trait/{gwas}".format(gwas=gwas)
	if os.path.exists("{}.results".format(fileout_prefix_ldsc_h2)):
		print("GWAS={} | LDSC outout file exists: {}. Will skip this LDSC regression...".format(gwas, fileout_prefix_ldsc_h2))
		continue
	file_gwas_sumstats = dict_gwas[gwas]["filepath"]
	
	if "pop-prev" in dict_gwas[gwas]:
		liability_string = "--pop-prev {} --samp-prev {}".format(dict_gwas[gwas]["pop-prev"], dict_gwas[gwas]["samp-prev"])
	else:
		liability_string = ""

	cmd = """{PYTHON2_EXEC} {flag_unbuffered} {script} --h2 {file_gwas_sumstats} \
    --ref-ld-chr {PATH_LDSC_DATA_MAIN}/baseline_v1.1_thin_annot/baseline. \
    --frqfile-chr {PATH_LDSC_DATA_MAIN}/1000G_Phase3_frq/1000G.EUR.QC. \
    --w-ld-chr {PATH_LDSC_DATA_MAIN}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --overlap-annot \
    --thin-annot \
    {liability_string} \
    --print-cov --print-coefficients --print-delete-vals \
    --out {fileout_prefix_ldsc_h2}""".format(
		PYTHON2_EXEC=PYTHON2_EXEC,
		flag_unbuffered="-u" if FLAG_UNBUFFERED else "",
		script=PATH_LDSC_SCRIPT,
		file_gwas_sumstats=file_gwas_sumstats,
		PATH_LDSC_DATA_MAIN=PATH_LDSC_DATA_MAIN,
		liability_string=liability_string,
		fileout_prefix_ldsc_h2=fileout_prefix_ldsc_h2
		)
	### NOT sure why "--print-cov" is needed. It outputs ".cov" file.		
	list_cmds_ldsc_prim.append(cmd)



### Call scheduler
job_scheduler(list_cmds=list_cmds_ldsc_prim, n_parallel_jobs=N_PARALLEL_LDSC_REGRESSION_JOBS)
# print(list_cmds_ldsc_prim)

###################################### OUTPUT --h2 ######################################

### Example
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.cov
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.delete
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.log
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.part_delete
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.results

print("Script is done!")



