
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import re
import sys

import string
import random



###################################### UTILS - get  ######################################


def get_cond_ref_ld_chr_name(cond_annotation, dataset):
	""" Function to get the ref_ld_chr_name for 'conditional annotation' for ldsc.py --h2/--h2-cts command """
	# *IMPORTANT*: ref_ld_chr_name MUST be full file path PLUS trailing "."
	dict_ldscore_path_prefix = {"mousebrain":"/scratch/sc-ldsc/celltypes.mousebrain.all/per_annotation/celltypes.mousebrain.all__mousebrain_all.{}.sem_mean.", # placeholder for .format() call works
						 		"tabula_muris":"/scratch/sc-ldsc/celltypes.tabula_muris.all/per_annotation/celltypes.tabula_muris.all__tabula_muris.{}.sem_mean.",
						 		}
	if not dataset in dict_ldscore_path_prefix:
		raise KeyError("dataset={} is not found in dict_ldscore_path_prefix.".format(dataset))
	cond_ref_ld_chr_name = dict_ldscore_path_prefix[dataset].format(cond_annotation)
	files_ldscore = glob.glob("{}*l2.ldscore.gz".format(cond_ref_ld_chr_name)) # get ldscore files for all chromosomes. glob() returns full file paths.
	if not len(files_ldscore) == 22: # we must have ldscore files for every chromosome, so the length 
		raise ValueError("dataset={} | cond_ref_ld_chr_name={} has n={} matching *l2.ldscore.gz files. Expected exactly 22 files. Check the ldscore file directory or update the dict_ldscore_path_prefix inside this function.".format(dataset, cond_ref_ld_chr_name, len(files_ldscore)))
	return(cond_ref_ld_chr_name)

# def get_list_cond_ref_ld_chr_name_LIST_VERSION(list_cond_annotations, dataset):
# 	""" Function to get the ref_ld_chr_name for 'conditional annotation' for ldsc.py --h2/--h2-cts command """
# 	# *IMPORTANT*: ref_ld_chr_name MUST be full file path PLUS trailing "."
# 	dict_ldscore_path_prefix = {"mousebrain":"/scratch/sc-ldsc/celltypes.mousebrain.all/per_annotation/celltypes.mousebrain.all__mousebrain_all.{}.sem_mean.", # placeholder for .format() call works
# 						 		"tabula_muris":"/scratch/sc-ldsc/celltypes.tabula_muris.all/per_annotation/celltypes.tabula_muris.all__tabula_muris.{}.sem_mean.",
# 						 		}
# 	if not dataset in dict_ldscore_path_prefix:
# 		raise KeyError("dataset={} is not found in dict_ldscore_path_prefix.".format(dataset))
# 	list_cond_ref_ld_chr_name = [dict_ldscore_path_prefix[dataset].format(x) for x in list_cond_annotations]
# 	for cond_ref_ld_chr_name in list_cond_ref_ld_chr_name:
# 		files_ldscore = glob.glob("{}*l2.ldscore.gz".format(cond_ref_ld_chr_name)) # get ldscore files for all chromosomes. glob() returns full file paths.
# 		if not len(files_ldscore) == 22: # we must have ldscore files for every chromosome, so the length 
# 			raise ValueError("dataset={} | cond_ref_ld_chr_name={} has n={} matching *l2.ldscore.gz files. Expected exactly 22 files. Check the ldscore file directory or update the dict_ldscore_path_prefix inside this function.".format(dataset, cond_ref_ld_chr_name, len(files_ldscore)))
# 	return(list_cond_ref_ld_chr_name)


###################################### UTILS - ALL GENES [COPY FROM MAIN] ######################################


def get_all_genes_ref_ld_chr_name(dataset):
	""" Function to get the ref_ld_chr_name for 'all genes annotation' for ldsc.py --h2/--h2-cts command """
	# *IMPORTANT*: ldsc_all_genes_ref_ld_chr_name MUST be full file path PLUS trailing "."
	dict_dataset_all_genes_path_prefix = {"mousebrain":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.mousebrain.",
						 				"tabula_muris":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.tabula_muris.",
						 				"campbell":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.campbell.",
						 				 }
	if not dataset in dict_dataset_all_genes_path_prefix:
		raise KeyError("dataset={} is not found in dict_dataset_all_genes_path_prefix.".format(dataset))
	ldsc_all_genes_ref_ld_chr_name = dict_dataset_all_genes_path_prefix[dataset]
	# some obnoxious validation of the matches
	files_ldscore = glob.glob("{}*l2.ldscore.gz".format(ldsc_all_genes_ref_ld_chr_name)) # get ldscore files for all chromosomes. glob() returns full file paths.
	if not len(files_ldscore) == 22: # we must have ldscore files for every chromosome, so the length 
		raise ValueError("dataset={} only has n={} matching {}*l2.ldscore.gz files. Expected 22 files. Check the ldscore file directory or update the dict_dataset_all_genes_path_prefix inside this function.".format(dataset, len(files_ldscore), ldsc_all_genes_ref_ld_chr_name))
	return(ldsc_all_genes_ref_ld_chr_name)


###################################### Job scheduler ######################################

def job_scheduler(list_cmds, n_parallel_jobs):
	""" Schedule parallel jobs with at most n_parallel_jobs parallel jobs."""
	list_of_processes = []
	batch = 1
	for i, cmd in enumerate(list_cmds, start=1):
		print("job schedule batch = {} | i = {} | Running command: {}".format(batch, i, cmd))
		## p = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
		### You need to keep devnull open for the entire life of the Popen object, not just its construction. 
		### FNULL = open(os.devnull, 'w') # devnull filehandle does not need to be closed?
		p = subprocess.Popen(cmd, shell=True)
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
###################################### PARAMS AND CONSTANTS ######################################
##################################################################################################



PYTHON3_EXEC = "/tools/anaconda/3-4.4.0/envs/py3_anaconda3_PT180510/bin/python3"
PYTHON2_EXEC = "/tools/anaconda/3-4.4.0/envs/py27_anaconda3_PT170705/bin/python2"

PATH_LDSC_SCRIPT = "/raid5/projects/timshel/sc-genetics/ldsc/ldsc-timshel/ldsc.py" 
N_PARALLEL_LDSC_REGRESSION_JOBS = 4

list_gwas = ["BMI_Yengo2018"]

# list_gwas = [
# "ADHD_PGC_Demontis2017",
# "AN_PGC_Duncan2017",
# "ASD_iPSYCH_PGC_Grove2018",
# "blood_EOSINOPHIL_COUNT",
# "EA3_Lee2018",
# "HEIGHT_Yengo2018",
# "LIPIDS_HDL_Willer2013",
# "MDD_PGC_Wray2018",
# "RA_Okada2014",
# "SCZ_Ripke2014",
# "WHR_adjBMI_Shungin2015",
# "WHR_Shungin2015",
# "INSOMNIA_Jansen2018",
# "BMI_Yengo2018",
# ]


################## Cell-types ##################

list_cond_annotations_mb = ["MEGLU10","MEGLU1","DEINH3","TEGLU23","MEINH2","MEGLU11","DEGLU5","TEGLU17"]
list_cond_annotations_tb = ["Brain_Non-Myeloid.neuron","Brain_Non-Myeloid.oligodendrocyte_precursor_cell"]

dict_cts_conditional = {"celltypes.mousebrain.all":
						{"file_cts":"/raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/celltypes.mousebrain.all.ldcts.txt",
						"dataset":"mousebrain",
						"conditional_annotations":list_cond_annotations_mb},
 					 "celltypes.tabula_muris.all":
 					  	{"file_cts":"/raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/celltypes.tabula_muris.all.ldcts.txt",
 					  	"dataset":"tabula_muris",
 					  	"conditional_annotations":list_cond_annotations_tb},
 					 }


#########################################################################################
###################################### PRE-PROCESS ######################################
#########################################################################################

#########################################################################################
###################################### RUN LDSC PRIM ######################################
#########################################################################################

### Create job commands
list_cmds_ldsc_prim = []
for prefix_genomic_annot, param_dict in dict_cts_conditional.items():
	ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(param_dict["dataset"])
	# list_cond_ref_ld_chr_name = get_list_cond_ref_ld_chr_name(param_dict["conditional_annotations"], param_dict["dataset"])
	for gwas in list_gwas:
		for cond_annotation in param_dict["conditional_annotations"]:
			cond_ref_ld_chr_name = get_cond_ref_ld_chr_name(cond_annotation, param_dict["dataset"])
			fileout_prefix = "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/{prefix_genomic_annot}__{gwas}__CONDITIONAL__{cond_annotation}".format(gwas=gwas, prefix_genomic_annot=prefix_genomic_annot, cond_annotation=cond_annotation)
			if os.path.exists("{}.cell_type_results.txt".format(fileout_prefix)):
				print("GWAS={}, prefix_genomic_annot={} | LDSC outout file exists: {}. Will skip this LDSC regression...".format(gwas, prefix_genomic_annot, fileout_prefix))
				continue
			### I'm 90% sure that ldsc.py ONLY runs on python2 - and not python3.
			### *OBS*: we are runnin ldsc python script with UNBUFFERED stdout and stderr
			### REF: https://stackoverflow.com/questions/230751/how-to-flush-output-of-print-function
			### python -u: Force the stdout and stderr streams to be unbuffered. THIS OPTION HAS NO EFFECT ON THE STDIN STREAM [or writing of other files, e.g. the ldsc .log file]. See also PYTHONUNBUFFERED.
			cmd = """{PYTHON2_EXEC} -u {script} --h2-cts /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/{gwas}.sumstats.gz \
			--ref-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/baseline_v1.1/baseline.,{ldsc_all_genes_ref_ld_chr_name},{cond_ref_ld_chr_name} \
			--w-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
			--ref-ld-chr-cts /raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/{prefix_genomic_annot}.ldcts.txt \
			--out {fileout_prefix}""".format(
				PYTHON2_EXEC=PYTHON2_EXEC,
				script=PATH_LDSC_SCRIPT,
				gwas=gwas,
				prefix_genomic_annot=prefix_genomic_annot,
				ldsc_all_genes_ref_ld_chr_name=ldsc_all_genes_ref_ld_chr_name,
				cond_ref_ld_chr_name=cond_ref_ld_chr_name,
				fileout_prefix=fileout_prefix
				)
			list_cmds_ldsc_prim.append(cmd)


### Call scheduler
job_scheduler(list_cmds=list_cmds_ldsc_prim, n_parallel_jobs=N_PARALLEL_LDSC_REGRESSION_JOBS)


###################################### XXXX ######################################


print("Script is done!")



