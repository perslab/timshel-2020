
from __future__ import print_function
import argparse
import glob
import os
import subprocess
import re

###################################### USAGE ######################################
# Compatibility: Python 2 and 3

###################################### TODO ######################################

### ARGS
# gwas
# prefix_genomic_annot file
# Outdir
# Flag to run without BASELINE
# Flag to run with/without all_genes --> give all genes file prefix


# ### Run the regression
# gwas=BMI_Yengo2018
# prefix_genomic_annot=novo_bulk_sema_lira
# python2 /raid5/projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py \
#     --h2-cts /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/${gwas}.sumstats.gz \
#     --ref-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline. \
#     --out /raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/wgcna_modules/${prefix_genomic_annot}_${gwas} \
#     --ref-ld-chr-cts /XXX/${prefix_genomic_annot}.ldcts \
#     --w-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights.

###################################### DESCRIPTION ######################################


### Output
# This script will call ldsc.py to do prefix_genomic_annot regression.
# The following output files will be written to the --prefix_annot_files:
# <OUT>.cell_type_results.txt
# <OUT>.log


###################################### WIKI ######################################

### DOCS weights and baseline
# --ref-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline.
# --ref-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/baseline_v1.1/baseline.

# --w-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights.
# --w-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.


###################################### FUNCTIONS ######################################
def preprocess(prefix_genomic_annot, file_multi_gene_set):
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
	cmd="{PYTHON_EXEC} wrapper_compute_ldscores.py --prefix_annot_files /scratch/sc-ldsc/{prefix_genomic_annot}/ --n_parallel_jobs 4".format(PYTHON_EXEC=PYTHON_EXEC, prefix_genomic_annot=prefix_genomic_annot)
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

	### make cts file
	###  *RESOURCE NOTES*: this script is light-weight and uses no computational resources
	cmd="{PYTHON_EXEC} make_cts_file.py --prefix_ldscore_files /scratch/sc-ldsc/{prefix_genomic_annot}/per_annotation/ --cts_outfile /raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/{prefix_genomic_annot}.ldcts.txt".format(PYTHON_EXEC=PYTHON_EXEC, prefix_genomic_annot=prefix_genomic_annot)
	# ^*OBS***:DIRTY USING  as prefix in  {prefix_genomic_annot}.ldcts.txt. FIX THIS.
	print("Running command: {}".format(cmd))
	p = subprocess.Popen(cmd, shell=True)
	p.wait()
	print("Return code: {}".format(p.returncode))
	# RUNTIME ----> 0 min
	if not p.returncode == 0:
		raise Exception("Got non zero return code")

######################################  ######################################

def get_all_genes_ref_ld_chr_name(prefix_genomic_annot, annot_name_all_genes="all_genes_in_dataset"):
	""" Function to get the ref_ld_chr_name for ldsc.py --h2/--h2-cts command """
	dir_per_annot = "/scratch/sc-ldsc/{prefix_genomic_annot}/per_annotation/".format(prefix_genomic_annot=prefix_genomic_annot) # *OBS* hardcoded
	files_per_annot = glob.glob("{}/*{}*".format(dir_per_annot,annot_name_all_genes)) # shortlist number files by globbing. Not needed, but more efficient
	if len(files_per_annot)==0:
		raise Exception("No files matching '{}' in dir_per_annot={}".format(annot_name_all_genes, dir_per_annot))
	dict_matches = {}
	for file_path in files_per_annot:
		m = re.search(r"^(.*%s.*)\.(\d{1,2})\.l2.ldscore.gz$" % annot_name_all_genes, os.path.basename(file_path)) # REF 'using a variable inside a regex' https://stackoverflow.com/a/6931048/6639640
		if m:
			ldsc_ref_ld_chr_name = m.groups()[0] # groups()[0]= celltypes.mousebrain.all.mousebrain_all.all_genes_in_dataset.dummy" if file_path= celltypes.mousebrain.all.mousebrain_all.all_genes_in_dataset.dummy.5.l2.ldscore.gz"
			chromosome = m.groups()[1]
			dict_matches[chromosome] = ldsc_ref_ld_chr_name
	# some obnoxious validation of the matches
	assert(set(dict_matches.keys()) == set(map(str, range(1,23)))) # we must have ldscore files for every chromosome
	assert(len(set(dict_matches.values())) == 1) # we must have only one unique 'basename' for ldsc_ref_ld_chr_name
	return os.path.join(dir_per_annot, ldsc_ref_ld_chr_name+".") # return full file path PLUS trailing "." which IS NEEDED

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



PYTHON_EXEC = "/tools/anaconda/2-4.4.0/bin/python2" # runs on python2
PATH_LDSC_SCRIPT = "/raid5/projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py" 
N_PARALLEL_LDSC_REGRESSION_JOBS = 3
FLAG_BINARY = True

list_gwas = ["BMI_Yengo2018"]
# list_gwas = ["BMI_Yengo2018",
# "EA3_Lee2018",
# "SCZ_Ripke2014",
# "HEIGHT_Yengo2018",
# "1KG_phase3_EUR_null_gwas_P1"]




################## Cell-types ##################
FLAG_WGCNA = False
# dict_genomic_annot = {"celltypes.mousebrain.all.mean":"/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/multi_geneset.mousebrain_all.mean.txt"}
dict_genomic_annot = {"celltypes.mousebrain.all.binary":"/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/multi_geneset.mousebrain_all.txt",
					  "celltypes.tabula_muris.all.binary":"/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/multi_geneset.tabula_muris.txt"}

# dict_genomic_annot = {"celltypes.campbell_lvl1.all":"/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/multi_geneset.campbell_lvl1.txt",
#					   "celltypes.campbell_lvl2.all":"/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/multi_geneset.campbell_lvl2.txt"}


################## WGCNA ##################
# FLAG_WGCNA = True
### Tabula muris
# dict_genomic_annot = {"wgcna.tabula_muris-181129":"/projects/jonatan/tmp-maca/tables/tabula_muris_1_cell_cluster_module_genes.csv"}


#########################################################################################
###################################### PRE-PROCESS ######################################
#########################################################################################

for prefix_genomic_annot, file_multi_gene_set in dict_genomic_annot.items():
	preprocess(prefix_genomic_annot, file_multi_gene_set)

#########################################################################################
###################################### RUN LDSC PRIM ######################################
#########################################################################################

### Create job commands
list_cmds_ldsc_prim = []
for prefix_genomic_annot in dict_genomic_annot:
	ldsc_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(prefix_genomic_annot)
	for gwas in list_gwas:
		fileout_prefix = "/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/{prefix_genomic_annot}.{gwas}.baseline_v1.1_all_genes".format(gwas=gwas, prefix_genomic_annot=prefix_genomic_annot)
		if os.path.exists("{}.cell_type_results.txt".format(fileout_prefix)):
			print("GWAS={}, prefix_genomic_annot={} | LDSC outout file exists: {}. Will skip this LDSC regression...".format(gwas, prefix_genomic_annot, fileout_prefix))
			continue
		cmd = """{PYTHON_EXEC} {script} --h2-cts /raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/{gwas}.sumstats.gz \
		--ref-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/baseline_v1.1/baseline.,{ldsc_ref_ld_chr_name} \
		--w-ld-chr /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
		--ref-ld-chr-cts /raid5/projects/timshel/sc-genetics/sc-genetics/src/ldsc/cts_files/{prefix_genomic_annot}.ldcts.txt \
		--out {fileout_prefix}""".format(
			PYTHON_EXEC=PYTHON_EXEC,
			script=PATH_LDSC_SCRIPT,
			gwas=gwas,
			prefix_genomic_annot=prefix_genomic_annot,
			ldsc_ref_ld_chr_name=ldsc_ref_ld_chr_name,
			fileout_prefix=fileout_prefix
			)
		list_cmds_ldsc_prim.append(cmd)


### Call scheduler
job_scheduler(list_cmds=list_cmds_ldsc_prim, n_parallel_jobs=N_PARALLEL_LDSC_REGRESSION_JOBS)


###################################### XXXX ######################################


print("Script is done!")



