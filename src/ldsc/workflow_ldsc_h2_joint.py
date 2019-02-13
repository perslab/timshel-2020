
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



###################################### UTILS - ALL GENES ######################################



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


list_gwas = ["BMI_UKBB_Loh2018"]

# list_gwas = ["BMI_UKBB_Loh2018",
# 			 "BMI_UPDATE_Yengo2018",
# 			 "T2D_DIAMANTE_Mahajan2018",
# 			 "T2D_UKBB_DIAMANTE_Mahajan2018",
# 			 "T2D_UKBB_Loh2018",
# 			 "HEIGHT_UKBB_Loh2018",
# 			 "HEIGHT_Yengo2018",
# 			 "LIPIDS_LDL_Teslovich2010",
# 			 "RA_Okada2014",
# 			]


##################################################################################################
############################################ PARAMS ##############################################
##################################################################################################

# ################## MOUSEBRAIN ##################
# dict_annotations = collections.defaultdict(dict)
# list_annotations = ["TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12"]
# for annotation in list_annotations:
# 	dict_annotations[annotation]["name_context"] = "mousebrain_all.{}.sem_mean".format(annotation)
# 	dict_annotations[annotation]["file_path_prefix"] = "/scratch/sc-ldsc/celltypes.mousebrain.all/per_annotation/celltypes.mousebrain.all__{}".format(dict_annotations[annotation]["name_context"])
# # ALT using dict comprehension: dict_annotations = {key:{"name_context":"mousebrain_all.{}.sem_mean".format(key)} for key in list_annotations}

# dict_annotations_mb = dict_annotations


# ################## TABULA MURIS ##################
# dict_annotations = collections.defaultdict(dict)
# list_annotations = ["Brain_Non-Myeloid.neuron","Brain_Non-Myeloid.oligodendrocyte_precursor_cell"]

# for annotation in list_annotations:
# 	dict_annotations[annotation]["name_context"] = "tabula_muris.{}.sem_mean".format(annotation)
# 	dict_annotations[annotation]["file_path_prefix"] = "/scratch/sc-ldsc/celltypes.tabula_muris.all/per_annotation/celltypes.tabula_muris.all__{}".format(dict_annotations[annotation]["name_context"])
# dict_annotations_tm = dict_annotations


# ################## COMBINE ##################
# # Output name: <runname/key_prim>__<annotation_id>__<gwas>
# dict_run = {"celltypes.mousebrain":
# 						{"dataset":"mousebrain",
# 						"dict_annotations":dict_annotations_mb},
# 			"celltypes.tabula_muris":
# 						{"dataset":"tabula_muris",
# 						"dict_annotations":dict_annotations_tm},
# 			}

# print(dict_run)


################## Cell-types ##################

# list_cond_annotations_mb = ["MEGLU10","MEGLU1","DEINH3","TEGLU23","MEINH2","MEGLU11","DEGLU5","TEGLU17"]
list_cond_annotations_mb = ["TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12"] # BMI_UKBB_Loh2018 FDR sign.
list_cond_annotations_tm = ["Brain_Non-Myeloid.neuron","Brain_Non-Myeloid.oligodendrocyte_precursor_cell"]


# format: annotation_id:cond_ref_ld_chr_name
# dict_addtional_cond_annotations_mb = {"wgcna.mousebrain-190111-dodgerblue":"/scratch/sc-ldsc/wgcna.mousebrain-190111.fdr_sign_celltypes.continuous/per_annotation/wgcna.mousebrain-190111.fdr_sign_celltypes.continuous__dodgerblue.", # must include trailing dot (".")} 

dict_cts_conditional = {"celltypes.mousebrain.all":
						{"ref_ld_chr_name":"DUMMY_TO_BE_UPDATED",
						"dataset":"mousebrain",
						"conditional_annotations":list_cond_annotations_mb},
 					 "celltypes.tabula_muris.all":
 					  	{"ref_ld_chr_name":"DUMMY_TO_BE_UPDATED",
 					  	"dataset":"tabula_muris",
 					  	"conditional_annotations":list_cond_annotations_tm},
 					 }


#########################################################################################
###################################### PRE-PROCESS ######################################
#########################################################################################


### add ref_ld_chr_name to dict
for run_name in list(dict_cts_conditional): # to avoid "RuntimeError: dictionary changed size during iteration”
	ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(dict_cts_conditional[run_name]["dataset"]) # just for testing that dataset exists....

	list_cond_ref_ld_chr_name = [get_cond_ref_ld_chr_name(cond_annotation, dict_cts_conditional[run_name]["dataset"]) for cond_annotation in dict_cts_conditional[run_name]["conditional_annotations"]] # get list of strings
	cond_ref_ld_chr_name = ",".join(list_cond_ref_ld_chr_name) # join by comma to parse this directly to LDSC
	annotation_id = "--".join(dict_cts_conditional[run_name]["conditional_annotations"])
	dict_cts_conditional[run_name]["ref_ld_chr_name"] = cond_ref_ld_chr_name

print(dict_cts_conditional)

#########################################################################################
###################################### RUN LDSC h2 ######################################
#########################################################################################


### About maximum size for a command single argument
# os.sysconf(os.sysconf_names["SC_ARG_MAX"])
# You can pass multiple strings of length 131071, but a single string arg cannot be longer than 131071 bytes.
# REF: https://stackoverflow.com/questions/29801975/why-is-the-subprocess-popen-argument-length-limit-smaller-than-what-the-os-repor

# /*
#  * These are the maximum length and maximum number of strings passed to the
#  * execve() system call.  MAX_ARG_STRLEN is essentially random but serves to
#  * prevent the kernel from being unduly impacted by misaddressed pointers.
#  * MAX_ARG_STRINGS is chosen to fit in a signed 32-bit integer.
#  */
# #define MAX_ARG_STRLEN (PAGE_SIZE * 32)
# #define MAX_ARG_STRINGS 0x7FFFFFFF

# ~$ echo $(( $(getconf PAGE_SIZE)*32 )) 
# 131072



### Create job commands
list_cmds_ldsc_prim = []
for run_name, param_dict in dict_cts_conditional.items():
	ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(param_dict["dataset"])
	annotation_id = "--".join(param_dict["conditional_annotations"]) # set name for output file. Use "--" as separator between annotation IDs. ALT "_-_" 
	cond_ref_ld_chr_name = param_dict["ref_ld_chr_name"]
	for gwas in list_gwas:
		# Output name: <runname/key_prim>__<annotation_id>__<gwas>
		fileout_prefix_ldsc_h2 = "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2_joint/{run_name}__{gwas}__JOINT__{annotation_id}".format(run_name=run_name, annotation_id=annotation_id, gwas=gwas)
		if os.path.exists("{}.results".format(fileout_prefix_ldsc_h2)):
			print("GWAS={}, run_name={},  annotation_id={} | LDSC outout file exists: {}. Will skip this LDSC regression...".format(gwas, run_name,  annotation_id, fileout_prefix_ldsc_h2))
			continue
		### I'm 90% sure that ldsc.py ONLY runs on python2 - and not python3.
		### *OBS*: we are runnin ldsc python script with UNBUFFERED stdout and stderr
		### REF: https://stackoverflow.com/questions/230751/how-to-flush-output-of-print-function
		### python -u: Force the stdout and stderr streams to be unbuffered. THIS OPTION HAS NO EFFECT ON THE STDIN STREAM [or writing of other files, e.g. the ldsc .log file]. See also PYTHONUNBUFFERED.
		cmd = """{PYTHON2_EXEC} {flag_unbuffered} {script} --h2 /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/{gwas}.sumstats.gz \
        --ref-ld-chr {PATH_LDSC_DATA_MAIN}/baseline_v1.1_thin_annot/baseline.,{ldsc_all_genes_ref_ld_chr_name},{cond_ref_ld_chr_name} \
        --frqfile-chr {PATH_LDSC_DATA_MAIN}/1000G_Phase3_frq/1000G.EUR.QC. \
        --w-ld-chr {PATH_LDSC_DATA_MAIN}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --overlap-annot \
        --thin-annot \
        --print-cov --print-coefficients --print-delete-vals \
        --out {fileout_prefix_ldsc_h2}""".format(
			PYTHON2_EXEC=PYTHON2_EXEC,
			flag_unbuffered="-u" if FLAG_UNBUFFERED else "",
			script=PATH_LDSC_SCRIPT,
			gwas=gwas,
			run_name=run_name,
			ldsc_all_genes_ref_ld_chr_name=ldsc_all_genes_ref_ld_chr_name,
			PATH_LDSC_DATA_MAIN=PATH_LDSC_DATA_MAIN,
			cond_ref_ld_chr_name=cond_ref_ld_chr_name,  # this is a file prefix. LDSC will add .<CHR>.l2.ldsc.gz to the file
			fileout_prefix_ldsc_h2=fileout_prefix_ldsc_h2
			)
		### NOT sure why "--print-cov" is needed. It outputs ".cov" file.

		### Ensure that command is not too long
		if len(cmd) >= 131072:
			print("Command length (number of characters): {}".format(len(cmd)))
			raise ValueError("Command length too long for linux to handle: cmd = {}".format(cmd))
		list_cmds_ldsc_prim.append(cmd)


### Call scheduler
job_scheduler(list_cmds=list_cmds_ldsc_prim, n_parallel_jobs=N_PARALLEL_LDSC_REGRESSION_JOBS)

###################################### OUTPUT --h2 ######################################

### Example
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.cov
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.delete
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.log
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.part_delete
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.results

print("Script is done!")



