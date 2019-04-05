
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import re
import sys

import string
import random

import make_annot_from_geneset_all_chr

###################################### TODO ######################################

### Run without the default use of two-stage estimation and a cap on maximum chisq in estimation of the intercept term
# --chisq-max 9999 --two-step 9999 # (an arbitrarily high value) per recommendation of the Neale Lab. This is meant to adjust for the very large sample size of the UK Biobank. 
# REF: http://www.nealelab.is/blog/2017/9/20/insights-from-estimates-of-snp-heritability-for-2000-traits-and-disorders-in-uk-biobank

###################################### USAGE ######################################
# Compatibility: Python 2 and 3

### Run in unbuffered mode
# time python -u workflow_ldsc_cts.py |& tee workflow_ldsc_cts.UNNAMED.out.txt

###################################### DOCS ######################################


###################################### DESCRIPTION ######################################


### Output
# This script will call ldsc.py to do prefix_genomic_annot regression.
# The following output files will be written to the --prefix_annot_files:
# <OUT>.cell_type_results.txt
# <OUT>.log


###################################### WIKI ######################################

### DOCS weights and baseline
# --ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline.
# --ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/baseline_v1.1/baseline.

# --w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights.
# --w-ld-chr /projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.


###################################### FUNCTIONS ######################################
def write_cts_file_filter(prefix_genomic_annot, file_multi_gene_set):
	""" 
	Write 'annotation filter' file based on annotations in file_multi_gene_set to use as input to make_cts_file.py. 
	This is to modify the default behavior of make_cts_file.py, which is to write a cts file for all annotations in /per_annotation/ directory.
	
	DESIGN [OUTDATED]:
	1) read log.{prefix_genomic_annot}.multi_geneset.txt file (generated by make_annot_from_geneset_all_chr.py)
		- this file ALWAYS contain the 3 columns: "annotation", "gene_input", "annotation_value"
		- this file ALWAYS contain ALL annotations in the out_dir.
	2) use unique annotations in file_multi_gene_set to filter file_multi_geneset_all.
	3) write file_cts_annotation_filter to /tmp/ dir 
	4) return file_cts_annotation_filter filename.
	"""
	print("Generating file_cts_annotation_filter...")
	str_random = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)) # geneate random ascii string of len 10. REF: https://stackoverflow.com/a/2257449/6639640
	file_cts_annotation_filter = "/tmp/{}.{}.txt".format(prefix_genomic_annot, str_random) # *OBS*: we add a random string to make sure any parallel processes are not trying to write to the same file as the same time.
	
	### Read existing multi_geneset containing ALL annotations
	# file_multi_geneset_all = "/scratch/sc-ldsc/{prefix_genomic_annot}/log.{prefix_genomic_annot}.multi_geneset.txt".format(prefix_genomic_annot=prefix_genomic_annot)
	# df_multi_geneset_all = pd.read_csv(file_multi_geneset_all, sep="\t", index_col=False)
	
	### Parse file_multi_gene_set
	df_multi_gene_set = make_annot_from_geneset_all_chr.read_multi_gene_set_file(file_multi_gene_set=file_multi_gene_set,
															 out_dir="/scratch/sc-ldsc/{prefix_genomic_annot}".format(prefix_genomic_annot=prefix_genomic_annot),
															 out_prefix=prefix_genomic_annot,
															 flag_encode_as_binary_annotation=False, # can be set to both True and False
															 flag_mouse=True if FLAG_WGCNA else False, 
															 flag_wgcna=True if FLAG_WGCNA else False,
															 print_log_files=False) # OBS: print_log_files=False is important
	annotations = df_multi_gene_set["annotation"].unique() # returns NumPy array.
	### Write file
	with open(file_cts_annotation_filter, "w") as fh_out:
		for annotation in annotations: 
			fh_out.write(annotation+"\n")
	print("Wrote n={} annotations to file_cts_annotation_filter={}.".format(len(annotations), file_cts_annotation_filter))
	return file_cts_annotation_filter


def ldsc_pre_computation(prefix_genomic_annot, file_multi_gene_set):

	### Error checks
	if "__" in prefix_genomic_annot:
		raise ValueError("prefix_genomic_annot={} contains double underscore ('__'). These characters are reserved keywords for splitting files downstream in the pipeline.".format(prefix_genomic_annot))

	### Make annot
	###  *RESOURCE NOTES*: if you have many modules (~3000-6000) then set --n_parallel_jobs to ~2-5 (instead of 22). Otherwise the script will up all the MEMORY on yggdrasil and fail.
	cmd = """{PYTHON3_EXEC} {flag_unbuffered} make_annot_from_geneset_all_chr.py \
	--file_multi_gene_set {file_multi_gene_set} \
	--file_gene_coord /projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
	--windowsize 100000 \
	--bimfile_basename /projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
	{flag_binary} \
	{flag_wgcna} \
	--out_dir /scratch/sc-ldsc/{prefix_genomic_annot} \
	--out_prefix {prefix_genomic_annot}
	""".format(PYTHON3_EXEC=PYTHON3_EXEC,
		flag_unbuffered="-u" if FLAG_UNBUFFERED else "", 
		file_multi_gene_set=file_multi_gene_set, 
		prefix_genomic_annot=prefix_genomic_annot, 
		flag_wgcna="--flag_wgcna --flag_mouse" if FLAG_WGCNA else "",
		flag_binary="--flag_encode_as_binary_annotation" if FLAG_BINARY else "",
		) 
	# --windowsize 100000 \ ---> 100 kb default
	# --n_parallel_jobs 11
	
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True, bufsize=0 if FLAG_UNBUFFERED else -1)
	p.wait()
	print("Return code: {}".format(p.returncode))
	if not p.returncode == 0:
		raise Exception("make_annot_from_geneset_all_chr.py: Got non zero return code running command:\n{}".format(cmd))


	### compute LD scores
	### *RESOURCE NOTES*: this script uses a lot of CPU. Never run more than 4 parallel jobs. 4 parallel jobs will use ~220% CPU
	cmd="{PYTHON3_EXEC} {flag_unbuffered} wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/{prefix_genomic_annot}/ --n_parallel_jobs 2".format(PYTHON3_EXEC=PYTHON3_EXEC,
																																								flag_unbuffered="-u" if FLAG_UNBUFFERED else "", 
																																								prefix_genomic_annot=prefix_genomic_annot)
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True, bufsize=0 if FLAG_UNBUFFERED else -1)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> ~6 h for ~500 modules with --n_parallel_jobs=4
	if not p.returncode == 0:
		raise Exception("wrapper_compute_ldscores.py: Got non zero return code running command:\n{}".format(cmd))

	### split LD scores
	### This script will read 1 ".COMBINED_ANNOT.$CHR.l2.ldscore.gz" file  (N_SNPs x N_ANNOTATION) per parallel process.
	###  *RESOURCE NOTES*: this script does not use much memory (it uses < 10-50GB?) and can easy be run with full parallelization (n=22)
	cmd="{PYTHON3_EXEC} {flag_unbuffered} split_ldscores.py --prefix_ldscore_files /scratch/sc-ldsc/{prefix_genomic_annot}/ --n_parallel_jobs 22".format(PYTHON3_EXEC=PYTHON3_EXEC,
																																								flag_unbuffered="-u" if FLAG_UNBUFFERED else "", 
																																								prefix_genomic_annot=prefix_genomic_annot)
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True, bufsize=0 if FLAG_UNBUFFERED else -1)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> ~10 min
	if not p.returncode == 0:
		raise Exception("split_ldscores.py: Got non zero return code running command:\n{}".format(cmd))

	### Write CTS filter
	file_cts_annotation_filter = write_cts_file_filter(prefix_genomic_annot, file_multi_gene_set)
	sys.stdout.flush()

	### make cts file
	###  *RESOURCE NOTES*: this script is light-weight and uses no computational resources
	cmd="{PYTHON3_EXEC} {flag_unbuffered} make_cts_file.py --prefix_ldscore_files /scratch/sc-ldsc/{prefix_genomic_annot}/per_annotation/ --cts_outfile /projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/{prefix_genomic_annot}.ldcts.txt --annotation_filter {file_cts_annotation_filter}".format(PYTHON3_EXEC=PYTHON3_EXEC,
																																																																								  flag_unbuffered="-u" if FLAG_UNBUFFERED else "", 
																																																																								  prefix_genomic_annot=prefix_genomic_annot,
																																																																								  file_cts_annotation_filter=file_cts_annotation_filter)
	# ^*OBS***:DIRTY USING  as prefix in  {prefix_genomic_annot}.ldcts.txt. FIX THIS.
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True, bufsize=0 if FLAG_UNBUFFERED else -1)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> 0 min
	if not p.returncode == 0:
		raise Exception("make_cts_file.py: Got non zero return code running command:\n{}".format(cmd))

###################################### UTILS - ALL GENES ######################################


def get_all_genes_ref_ld_chr_name(dataset):
	""" Function to get the ref_ld_chr_name for 'all genes annotation' for ldsc.py --h2/--h2-cts command """
	# *IMPORTANT*: ldsc_all_genes_ref_ld_chr_name MUST be full file path PLUS trailing "."
	dict_dataset_all_genes_path_prefix = {"mousebrain":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.mousebrain.",
						 				"tabula_muris":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.tabula_muris.",
						 				"campbell":"/scratch/sc-ldsc/control.all_genes_in_dataset/per_annotation/control.all_genes_in_dataset__all_genes_in_dataset.campbell.",
						 				 "dataset_with_no_all_genes":"" # value must be empty string.
						 				 }
	if not dataset in dict_dataset_all_genes_path_prefix:
		raise KeyError("dataset={} is not found in dict_dataset_all_genes_path_prefix.".format(dataset))
	ldsc_all_genes_ref_ld_chr_name = dict_dataset_all_genes_path_prefix[dataset]
	if ldsc_all_genes_ref_ld_chr_name: # only needed to support dataset_with_no_all_genes (empty string valuates false)
		# some obnoxious validation of the matches
		files_ldscore = glob.glob("{}*l2.ldscore.gz".format(ldsc_all_genes_ref_ld_chr_name)) # get ldscore files for all chromosomes. glob() returns full file paths.
		if not len(files_ldscore) == 22: # we must have ldscore files for every chromosome, so the length 
			raise ValueError("dataset={} only has n={} matching {}*l2.ldscore.gz files. Expected 22 files. Check the ldscore file directory or update the dict_dataset_all_genes_path_prefix inside this function.".format(dataset, len(files_ldscore), ldsc_all_genes_ref_ld_chr_name))
	return(ldsc_all_genes_ref_ld_chr_name)



# def get_all_genes_ref_ld_chr_name_V1(prefix_genomic_annot, annot_name_all_genes="all_genes_in_dataset"):
# 	""" Function to get the ref_ld_chr_name for 'all genes annotation' for ldsc.py --h2/--h2-cts command """
# 	dir_per_annot = "/scratch/sc-ldsc/{prefix_genomic_annot}/per_annotation/".format(prefix_genomic_annot=prefix_genomic_annot) # *OBS* hardcoded
# 	files_per_annot = glob.glob("{}/*{}*".format(dir_per_annot,annot_name_all_genes)) # shortlist number files by globbing. Not needed, but more efficient
# 	if len(files_per_annot)==0:
# 		raise Exception("No files matching '{}' in dir_per_annot={}".format(annot_name_all_genes, dir_per_annot))
# 	dict_matches = {}
# 	for file_path in files_per_annot:
# 		m = re.search(r"^(.*%s.*)\.(\d{1,2})\.l2.ldscore.gz$" % annot_name_all_genes, os.path.basename(file_path)) # REF 'using a variable inside a regex' https://stackoverflow.com/a/6931048/6639640
# 		if m:
# 			ldsc_all_genes_ref_ld_chr_name = m.groups()[0] # groups()[0]= celltypes.mousebrain.all.mousebrain_all.all_genes_in_dataset.dummy" if file_path= celltypes.mousebrain.all.mousebrain_all.all_genes_in_dataset.dummy.5.l2.ldscore.gz"
# 			chromosome = m.groups()[1]
# 			dict_matches[chromosome] = ldsc_all_genes_ref_ld_chr_name
# 	# some obnoxious validation of the matches
# 	assert(set(dict_matches.keys()) == set(map(str, range(1,23)))) # we must have ldscore files for every chromosome
# 	assert(len(set(dict_matches.values())) == 1) # we must have only one unique 'basename' for ldsc_all_genes_ref_ld_chr_name
# 	return os.path.join(dir_per_annot, ldsc_all_genes_ref_ld_chr_name+".") # return full file path PLUS trailing "." which IS NEEDED

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
###################################### PARAMS AND CONSTANTS ######################################
##################################################################################################



PYTHON3_EXEC = "/tools/anaconda/3-4.4.0/envs/py3_anaconda3_PT180510/bin/python3"
PYTHON2_EXEC = "/tools/anaconda/3-4.4.0/envs/py27_anaconda3_PT170705/bin/python2"

PATH_LDSC_SCRIPT = "/projects/timshel/sc-genetics/ldsc/ldsc-timshel/ldsc.py" 
FLAG_UNBUFFERED = True
N_PARALLEL_LDSC_REGRESSION_JOBS = 1
# FLAG_BINARY = True
FLAG_BINARY = False


# list_gwas = ["BMI_UKBB_Loh2018"] # BMI_UPDATE_Yengo2018

# list_gwas = ["BMI_UKBB_Loh2018_no_mhc_max_chisq_80",
# "BMI_UKBB_Loh2018_no_mhc_max_chisq_720",
# "BMI_UPDATE_Yengo2018_no_mhc_max_chisq_720",
# "BMI_UPDATE_Yengo2018_no_mhc_max_chisq_80"]


list_gwas = ["1KG_phase3_EUR_null_gwas_P{}".format(x) for x in range(1,11)] # 1..10

# list_gwas = ["ADHD_PGC_Demontis2017",
# "AD_Jansen2019",
# "AD_Lambert2013",
# "AN_PGC_Duncan2017",
# "ASD_iPSYCH_PGC_Grove2018",
# "BIP_PGC2018",
# "blood_EOSINOPHIL_COUNT",
# "BMI_Locke2015",
# "BMI_UPDATE_Yengo2018",
# "CAD_Schunkert2011",
# "CELIAC_Dubois2010",
# "CROHNS_Jostins2012",
# "DEPRESSED_AFFECT_Nagel2018",
# "DEPRESSION_Nagel2018",
# "DS_Okbay2016",
# "EA2_Okbay2016",
# "EA3_Lee2018",
# "FG_Female_Lagou2018",
# "FG_Male_Lagou2018",
# "FI_Female_Lagou2018",
# "FI_Male_Lagou2018",
# "HBA1C_MAGIC_Wheeler2017",
# "LIPIDS_HDL_Teslovich2010",
# "HEIGHT_Wood2014",
# "HEIGHT_Yengo2018",
# "IBD_Jostins2012",
# "INSOMNIA_Jansen2018",
# "INTELLIGENCE_Savage2018",
# "INTELLIGENCE_Sniekers2017",
# "LIPIDS_LDL_Teslovich2010",
# "LUPUS_2015",
# "MDD_PGC_Wray2018",
# "MS_Patsopoulos2011",
# "NEUROTICISM_Nagel2018",
# "NEUROTICISM_OKBAY2016",
# "PBC_Cordell2015",
# "RA_Okada2014",
# "RB_Linner_2019",
# "SCZ_Pardinas2018",
# "SCZ_EUR_Ripke2014",
# "SWB_Okbay2016",
# "T1D_Bradfield2011",
# "T2DadjBMI_DIAMANTE_Mahajan2018",
# "T2D_DIAMANTE_Mahajan2018",
# "T2D_UKBB_DIAMANTE_Mahajan2018",
# "T2D_Xue2018",
# "LIPIDS_TG_Teslovich2010",
# "UC_Jostins2012",
# "WHR_adjBMI_Shungin2015",
# "WHR_Shungin2015",
# "WORRY_Nagel2018",
# "WHR_Pulit2019",
# "WHRadjBMI_Pulit2019",
# "BMI_Pulit2019",
# "BMI_male_Pulit2019",
# "BMI_female_Pulit2019",
# "BMI_UKBB_Loh2018",
# "WHRadjBMI_UKBB_Loh2018",
# "HEIGHT_UKBB_Loh2018",
# "DIASTOLICadjMED_UKBB_Loh2018",
# "SYSTOLICadjMED_UKBB_Loh2018",
# "CARDIOVASCULAR_UKBB_Loh2018",
# "MDD_Howard2019",
# "LIPIDS_HDL_Willer2013",
# "LIPIDS_LDL_Willer2013",
# "LIPIDS_TG_Willer2013",
# "LIPIDS_TC_Willer2013",
# "T2D_UKBB_Loh2018"]





################## Cell-types ##################
FLAG_WGCNA = False


### Mean MB [mousebrain_all_190306_es_fix]
# dict_genomic_annot = {"celltypes.mousebrain_190306_es_fix.all":
# 						{"dataset":"mousebrain",
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.mousebrain_all_190306_es_fix.txt.gz"},
#  					 }


### Mousebrain hierarchical (11 FDR cell-types + neurons)
# dict_genomic_annot = {"celltypes.mousebrain.bmi_loh2018_11fdr_celltypes":
# 						{"dataset":"mousebrain",
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.mousebrain_bmi_loh2018_11fdr_celltypes.sem_mean.txt.gz"},
#  					 "celltypes.mousebrain.neurons":
#  					  	{"dataset":"mousebrain",
#  					  	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.mousebrain_neurons.sem_mean.txt.gz"}
#  					 }


# ### Adipocyte (sem_mean only)
# dict_genomic_annot = {"celltypes.preadipocyte_developing_1808_branch":
# 						{"dataset":"dataset_with_no_all_genes",
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.preadipocyte_developing_1808_branch.sem_mean.txt.gz"},
#  					 "celltypes.preadipocyte_developing_1808_branch_pc2_quantile":
#  					  	{"dataset":"dataset_with_no_all_genes",
#  					  	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.preadipocyte_developing_1808_branch_pc2_quantile.sem_mean.txt.gz"}
#  					 }

# ### Mean TM ONLY [for null]
dict_genomic_annot = {"celltypes.tabula_muris.all":
 					  	{"dataset":"tabula_muris",
 					  	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.tabula_muris.sem_mean.txt"}
 					 }

# # ### Mean MB+TM
# dict_genomic_annot = {"celltypes.mousebrain.all":
# 						{"dataset":"mousebrain",
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.mousebrain_all.sem_mean.txt"},
#  					 "celltypes.tabula_muris.all":
#  					  	{"dataset":"tabula_muris",
#  					  	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.tabula_muris.sem_mean.txt"}
#  					 }


### top10pct (Skene and Hillary)
# dict_genomic_annot = {"celltypes.mousebrain_top10pct.all":
# 						{"dataset":"mousebrain",
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.mousebrain_all_top10pct_binary.txt.gz"},
#  					 "celltypes.tabula_muris_top10pct.all":
#  					  	{"dataset":"tabula_muris",
#  					  	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.tabula_muris_top10pct_binary.txt.gz"}
#  					 }

# ### Raw SEMs
# dict_genomic_annot = {"celltypes.mousebrain_raw_sems.all":
# 						{"dataset":"mousebrain",
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.mousebrain_all_raw_sems.txt.gz"},
#  					 "celltypes.tabula_muris_raw_sems.all":
#  					  	{"dataset":"tabula_muris",
#  					  	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.tabula_muris_raw_sems.txt.gz"}
#  					 }


# dict_genomic_annot = {"celltypes.campbell_lvl1.all":
# 						{"dataset":"campbell",
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.campbell_lvl1.sem_mean.txt"},
# 					 "celltypes.campbell_lvl2.all":
# 					   	{"dataset":"campbell",
# 					   	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/multi_geneset.campbell_lvl2.sem_mean.txt"}
# 					  }



### CMD find . -name "multi_geneset.*.txt" | xargs -I {} sh -c "basename {} .txt" | xargs -I {} sh -c "egrep 'sem_mean|all_genes_in_dataset' {}.txt > {}.sem_mean.txt"
### CMD [simple, creates double .txt]: find . -name "multi_geneset.*.txt" | xargs -I {} sh -c "egrep 'sem_mean|all_genes_in_dataset' {} > {}.sem_mean.txt"


################## WGCNA ##################
# FLAG_WGCNA = True


### Modules from FDR significant cell-types [v4, 190218] | WGCNA med deepSplit 1 + n=11 BMI_UKBB_Loh2018 FDR cell-types 
# Feb 18 4:22 PM
# - koerslen hedder Neurons_sub_ClusterName_7.3_run1
# - data.frame med module genes: /projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.3_run1_cell_cluster_module_genes.csv.gz
# - kMEs: /projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.3_run1_kMs_full_join.csv.gz
# Nogle observationer:
# - Baade lavenderblush og lightpink3 er 100% bevarede (nu som henholdsvis ‘palevioletred’ og ‘brown’ - beklager, har ikke ville begynde at rette i det gamle script :S)
# - antal moduler er faldet til 126

# dict_genomic_annot = {"wgcna.mousebrain-190218.fdr_sign_celltypes.continuous": 
# 					  	{"dataset":"mousebrain", 
# 					  	"file_multi_gene_set":"/projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.3_run1_cell_cluster_module_genes.csv.gz"}
# 					 }


### Mousebrain - trait specificity for lavenderblush and lightpink3 [190218]
# dict_genomic_annot = {"wgcna.mousebrain-190213.fdr_sign_celltypes.continuous": 
# 					  	{"dataset":"mousebrain", 
# 					  	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/data/gene_lists/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes--lavenderblush--lightpink3.csv"}
# 					 }

### Modules from FDR significant cell-types [v3, 190213] | NEW: WGCNA run on n=11 BMI_UKBB_Loh2018 FDR cell-types 
# Feb 13 [6:15 PM]
# Hej Pascal,
# Nu er den nye mb WGCNA analyse  endelig klar. Jeg gik tilbage til den tidligere udgave af pipelinen fra Januar.
# Koerslen hedder Neurons_sub_ClusterName_7.2_run1 og den relevante output fil er
# /projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz
# Dodgerblue modulet gaar igen uden aendring, nu under navnet “lavenderblush”.
# (har i senere aendringer i scriptet soerget for at alle random seeds, inklusive navne, er reproducible) 
# dict_genomic_annot = {"wgcna.tabula_muris-190111.fdr_sign_celltypes.continuous": 
# 						{"dataset":"tabula_muris", 
# 						"file_multi_gene_set":"/projects/jonatan/tabula_muris_3/tables/tabula_muris_3_cell_cluster_module_genes.csv.gz"},
# 					  "wgcna.mousebrain-190213.fdr_sign_celltypes.continuous": 
# 					  	{"dataset":"mousebrain", 
# 					  	"file_multi_gene_set":"/projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz"}
# 					 }


# ### Modules from FDR significant cell-types [v2, 190111] + KME CUT-OFF
# ### NOTES kME reassign disabled; deepSplit = 2
# dict_genomic_annot = {}
# for cutoff_kme in [0.25, 0.30, 0.40, 0.50]:
# 	for dataset_name in ["tabula_muris", "mousebrain"]:
# 		key_name = "wgcna.{}-190111.fdr_sign_celltypes_min_kme_{:.0f}.continuous".format(dataset_name, cutoff_kme*100) # *100 to avoid dots in name
# 		if dataset_name == "tabula_muris": 
# 			file_prefix_cluster = "tabula_muris_3_cell_cluster_module_genes"
# 		elif dataset_name == "mousebrain": 
# 			file_prefix_cluster = "mb_Neurons_ClusterName_7_cell_cluster_module_genes"
# 		else:
# 			raise Exception(".")
# 		dict_genomic_annot[key_name] = {"dataset":dataset_name, 
# 									"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/data/gene_lists/{}.min_kme_{:.2f}.csv.gz".format(file_prefix_cluster, cutoff_kme)
# 									}


# ### Modules from FDR significant cell-types [v2, 190111]
# ### NOTES kME reassign disabled; deepSplit = 2
# dict_genomic_annot = {"wgcna.tabula_muris-190111.fdr_sign_celltypes.continuous": 
# 						{"dataset":"tabula_muris", 
# 						"file_multi_gene_set":"/projects/jonatan/tabula_muris_3/tables/tabula_muris_3_cell_cluster_module_genes.csv.gz"},
# 					  "wgcna.mousebrain-190111.fdr_sign_celltypes.continuous": 
# 					  	{"dataset":"mousebrain", 
# 					  	"file_multi_gene_set":"/projects/jonatan/mousebrain_7/tables/mb_Neurons_ClusterName_7_cell_cluster_module_genes.csv.gz"}
# 					 }



# ### Modules from FDR significant cell-types [v1, 181214]
# ### NOTES: too many modules generated. ~50 modules per cell-type
# dict_genomic_annot = {"wgcna.tabula_muris-181214.fdr_sign_celltypes.continuous": 
# 						{"dataset":"tabula_muris", 
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/data/gene_lists/tabula_muris-181214.tabula_muris_2_cell_cluster_module_genes.fdr_sign_celltypes.csv"},
# 					  "wgcna.mousebrain-181214.fdr_sign_celltypes.continuous": 
# 					  	{"dataset":"mousebrain", 
# 					  	"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/data/gene_lists/mousebrain-181214.mb_ClusterName_6_cell_cluster_module_genes.fdr_sign_celltypes.csv"}
# 					 }

### All modules
# dict_genomic_annot = {"wgcna.tabula_muris-181214.continuous":"/projects/jonatan/tabula_muris_2/tables/tabula_muris_2_cell_cluster_module_genes.csv.gz",
# 					  "wgcna.mousebrain-181214.continuous":"/projects/jonatan/tmp_mousebrain_6/tables/mb_ClusterName_6_cell_cluster_module_genes.csv.gz"}

### TEST
# dict_genomic_annot = {"wgcna.tmp_test_file_multi_gene_set_wgcna200": 
# 						{"dataset":"tabula_muris", 
# 						"file_multi_gene_set":"/projects/timshel/sc-genetics/sc-genetics/src/ldsc/multi_geneset_files/test_file_multi_gene_set_wgcna200.csv"},
# 					}



#########################################################################################
###################################### PRE-PROCESS ######################################
#########################################################################################

### Make sure that all_genes annotation is present, so program does not fail later on.
for prefix_genomic_annot in dict_genomic_annot.keys():
	param_dict = dict_genomic_annot[prefix_genomic_annot]
	ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(param_dict["dataset"])

### Run pre-computation
for prefix_genomic_annot in list(dict_genomic_annot.keys()): # list() needed for py3 compatibility. REF: https://stackoverflow.com/a/11941855/6639640
	param_dict = dict_genomic_annot[prefix_genomic_annot]
	try:
		ldsc_pre_computation(prefix_genomic_annot, param_dict["file_multi_gene_set"])
	except Exception as e:
		print("Caught exception during ldsc_pre_computation for prefix_genomic_annot={}.".format(prefix_genomic_annot))
		print("Exception: {}".format(e))
		print("Will drop prefix_genomic_annot={} from dict_genomic_annot and not do any further computations on this prefix_genomic_annot.".format(prefix_genomic_annot))
		dict_genomic_annot.pop(prefix_genomic_annot, None) # drop key from dict while iterating over it. REF: https://stackoverflow.com/questions/5384914/how-to-delete-items-from-a-dictionary-while-iterating-over-it and https://stackoverflow.com/a/11277439/6639640

#########################################################################################
###################################### RUN LDSC PRIM ######################################
#########################################################################################

### Create job commands
list_cmds_ldsc_prim = []
for prefix_genomic_annot, param_dict in dict_genomic_annot.items():
	ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(param_dict["dataset"])
	flag_all_genes = True
	if ldsc_all_genes_ref_ld_chr_name=="":
		print("OBS: Running without all_genes correction.")
		flag_all_genes = False
	for gwas in list_gwas:
		fileout_prefix = "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/{prefix_genomic_annot}__{gwas}".format(gwas=gwas, prefix_genomic_annot=prefix_genomic_annot)
		if os.path.exists("{}.cell_type_results.txt".format(fileout_prefix)):
			print("GWAS={}, prefix_genomic_annot={} | LDSC outout file exists: {}. Will skip this LDSC regression...".format(gwas, prefix_genomic_annot, fileout_prefix))
			continue
		### I'm 90% sure that ldsc.py ONLY runs on python2 - and not python3.
		### *OBS*: we are runnin ldsc python script with UNBUFFERED stdout and stderr
		### REF: https://stackoverflow.com/questions/230751/how-to-flush-output-of-print-function
		### python -u: Force the stdout and stderr streams to be unbuffered. THIS OPTION HAS NO EFFECT ON THE STDIN STREAM [or writing of other files, e.g. the ldsc .log file]. See also PYTHONUNBUFFERED.
		cmd = """{PYTHON2_EXEC} {flag_unbuffered} {script} --h2-cts /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/{gwas}.sumstats.gz \
		--ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/baseline_v1.1/baseline.{flag_all_genes}{ldsc_all_genes_ref_ld_chr_name} \
		--w-ld-chr /projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
		--ref-ld-chr-cts /projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/{prefix_genomic_annot}.ldcts.txt \
		--out {fileout_prefix}""".format(
			PYTHON2_EXEC=PYTHON2_EXEC,
			flag_unbuffered="-u" if FLAG_UNBUFFERED else "",
			script=PATH_LDSC_SCRIPT,
			gwas=gwas,
			prefix_genomic_annot=prefix_genomic_annot,
			flag_all_genes="," if flag_all_genes else "",
			ldsc_all_genes_ref_ld_chr_name=ldsc_all_genes_ref_ld_chr_name,
			fileout_prefix=fileout_prefix
			)
		list_cmds_ldsc_prim.append(cmd)


### Call scheduler
job_scheduler(list_cmds=list_cmds_ldsc_prim, n_parallel_jobs=N_PARALLEL_LDSC_REGRESSION_JOBS)


###################################### XXXX ######################################


print("Script is done!")



