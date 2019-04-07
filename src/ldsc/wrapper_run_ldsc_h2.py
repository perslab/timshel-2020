
from __future__ import print_function
import argparse
import glob
import os
import subprocess

###################################### USAGE ######################################
# Compatibility: Python 2 and 3

###################################### TODO ######################################

###################################### DESCRIPTION ######################################

# Run --h2 command (not --h2-cts). Estimating heritability enrichment (i.e., %h2/%SNPs) of any annotation.

###################################### CONSTANTS ######################################


python_exec = "/tools/anaconda/2-4.4.0/bin/python2" # runs on python2
path_ldsc_script = "/projects/timshel/sc-genetics/ldsc/ldsc/ldsc.py" 
# w-ld-chr
# ref-ld-chr --> baseline

###################################### ARGS ######################################

n_parallel_jobs = 20

annotations = ["Adipose_Subcutaneous",
"Adipose_Visceral_(Omentum)",
"Adrenal_Gland",
"Artery_Aorta",
"Artery_Coronary",
"Artery_Tibial",
"Bladder",
"Brain_Amygdala",
"Brain_Anterior_cingulate_cortex_(BA24)",
"Brain_Caudate_(basal_ganglia)",
"Brain_Cerebellar_Hemisphere",
"Brain_Cerebellum",
"Brain_Cortex",
"Brain_Frontal_Cortex_(BA9)",
"Brain_Hippocampus",
"Brain_Hypothalamus",
"Brain_Nucleus_accumbens_(basal_ganglia)",
"Brain_Putamen_(basal_ganglia)",
"Brain_Spinal_cord_(cervical_c-1)",
"Brain_Substantia_nigra",
"Breast_Mammary_Tissue",
"Cells_EBV-transformed_lymphocytes",
"Cells_Transformed_fibroblasts",
"Cervix_Ectocervix",
"Cervix_Endocervix",
"Colon_Sigmoid",
"Colon_Transverse",
"Esophagus_Gastroesophageal_Junction",
"Esophagus_Mucosa",
"Esophagus_Muscularis",
"Fallopian_Tube",
"Heart_Atrial_Appendage",
"Heart_Left_Ventricle",
"Kidney_Cortex",
"Liver",
"Lung",
"Minor_Salivary_Gland",
"Muscle_Skeletal",
"Nerve_Tibial",
"Ovary",
"Pancreas",
"Pituitary",
"Prostate",
"Skin_Not_Sun_Exposed_(Suprapubic)",
"Skin_Sun_Exposed_(Lower_leg)",
"Small_Intestine_Terminal_Ileum",
"Spleen",
"Stomach",
"Testis",
"Thyroid",
"Uterus",
"Vagina"] # this order it taken from the GTEX CTS files. MUST be kept in correct order. 
### "Index" = 1..52.
### Adipose_Subcutaneous = .../Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.1.
### Vagina = .../Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.52.


###################################### MAIN ######################################

# parser = argparse.ArgumentParser()
# parser.add_argument('--n_parallel_jobs', type=int, default=1, help='Number of parallel jobs to run.')
# args = parser.parse_args()


###################################### CALL ######################################

list_cmds = []
for annotation_index in range(1,len(annotations)+1): # 1...52, Hilary annotations are 1-based
	annotation_name = annotations[annotation_index-1] # python arrays are 0-based

	### -baseline, -all_genes
	cmd = """{python_exec} {script} \
	--h2 /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018.sumstats.gz \
	--ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.{annotation_index}. \
	--out '/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment-h2.no_baseline.no_all_genes/experiment.gtex_tissues.h2.{annotation_name}.BMI_Yengo2018' \
	--w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. \
	--frqfile-chr /projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_frq/1000G.EUR.QC. \
	--overlap-annot \
	--print-coefficients""".format(
		python_exec=python_exec,
		script=path_ldsc_script,
		annotation_index=annotation_index,
		annotation_name=annotation_name)
	list_cmds.append(cmd)

	### +baseline, -all_genes
	cmd = """{python_exec} {script} \
	--h2 /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018.sumstats.gz \
	--ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.{annotation_index}.,/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline. \
	--out '/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment-h2.baseline.no_all_genes/experiment.gtex_tissues.h2.{annotation_name}.BMI_Yengo2018' \
	--w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. \
	--frqfile-chr /projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_frq/1000G.EUR.QC. \
	--overlap-annot \
	--print-coefficients""".format(
		python_exec=python_exec,
		script=path_ldsc_script,
		annotation_index=annotation_index,
		annotation_name=annotation_name)
	list_cmds.append(cmd)

	### -baseline, +all_genes
	cmd = """{python_exec} {script} \
	--h2 /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018.sumstats.gz \
	--ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.{annotation_index}.,/projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control. \
	--out '/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment-h2.no_baseline.all_genes/experiment.gtex_tissues.h2.{annotation_name}.BMI_Yengo2018' \
	--w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. \
	--frqfile-chr /projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_frq/1000G.EUR.QC. \
	--overlap-annot \
	--print-coefficients""".format(
		python_exec=python_exec,
		script=path_ldsc_script,
		annotation_index=annotation_index,
		annotation_name=annotation_name)
	list_cmds.append(cmd)
	
	### +baseline, +all_genes
	cmd = """{python_exec} {script} \
	--h2 /projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_ldsc/timshel-collection/BMI_Yengo2018.sumstats.gz \
	--ref-ld-chr /projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.{annotation_index}.,/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_baseline/baseline.,/projects/timshel/sc-genetics/ldsc/data/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores/GTEx.control. \
	--out '/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment-h2.baseline.all_genes/experiment.gtex_tissues.h2.{annotation_name}.BMI_Yengo2018' \
	--w-ld-chr /projects/timshel/sc-genetics/ldsc/data/weights_hm3_no_hla/weights. \
	--frqfile-chr /projects/timshel/sc-genetics/ldsc/data/1000G_Phase3_frq/1000G.EUR.QC. \
	--overlap-annot \
	--print-coefficients""".format(
		python_exec=python_exec,
		script=path_ldsc_script,
		annotation_index=annotation_index,
		annotation_name=annotation_name)
	list_cmds.append(cmd)




### You need to keep devnull open for the entire life of the Popen object, not just its construction. 
# FNULL = open(os.devnull, 'w') # devnull filehandle does not need to be closed?

list_of_processes = []
batch = 1
for i, cmd in enumerate(list_cmds, start=1):
	print("batch = {} | i = {} | Running command: {}".format(batch, i, cmd))
	# p = subprocess.Popen(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	p = subprocess.Popen(cmd, shell=True)
	list_of_processes.append(p)
	print("batch = {} | i = {} | PIDs of running jobs (list_of_processes):".format(batch, i))
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



print("Script is done!")



