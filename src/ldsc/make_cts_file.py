
from __future__ import print_function
import argparse
import glob
import os
import subprocess

import re
import collections

###################################### USAGE ######################################
# Compatibility: Python 2 and 3

# python make_cts_file.py --prefix_ldscore_files /scratch/sc-ldsc/nn_lira_sema/nn_lira_sema --cts_outfile wgcna.nn_lira_sema.ldcts.txt # ---> FAILTS
# python make_cts_file.py --prefix_ldscore_files /scratch/sc-ldsc/nn_lira_sema-combined/per_annotation/nn_lira_sema --cts_outfile wgcna.nn_lira_sema.ldcts.txt
# python make_cts_file.py --prefix_ldscore_files /scratch/sc-ldsc/maca/per_annotation/ --cts_outfile wgcna.maca.ldcts.txt
# python make_cts_file.py --prefix_ldscore_files /scratch/sc-ldsc/mousebrain/per_annotation/ --cts_outfile wgcna.mousebrain.ldcts.txt
# python make_cts_file.py --prefix_ldscore_files /scratch/sc-ldsc/hypothalamus_mette_thesis/per_annotation/ --cts_outfile wgcna.hypothalamus_mette_thesis.ldcts.txt


###################################### TODO ######################################

# MAYBE Write CTS file to outdir (to keep it as a log)
# 

### Args to add
# Input dir
# Pattern (comma seperated) to match filename for which annotations to run on
# 


###################################### DESCRIPTION ######################################

### Purpose
# Make CTS file

### Steps

### Output



###################################### MAIN ######################################

parser = argparse.ArgumentParser()
parser.add_argument('--prefix_ldscore_files', type=str, help="""File path prefix to LDSC ldscore files to write out CTS file for. 
	Will make cts files for ldscores files matching <prefix_ldscore_files>*.l2.ldscore.gz
	Notice that if prefix_ldscore_files is a directory (with a trailing '/'), then all ldscore file in the directory will be used.
	Example with prefix as input: /scratch/sc-ldsc/ldsc_analysis_dir/files_prefix
	Example with dir as input: /scratch/sc-ldsc/ldsc_analysis_dir/ [OBS the trailing slash is important!]""")
parser.add_argument('--cts_outfile', type=str, help="Filename for output CTS file.")


args = parser.parse_args()


file_ldscores = glob.glob("{}*.l2.ldscore.gz".format(args.prefix_ldscore_files)) # e.g. <SOMEPATH>/nn_lira_sema.yellow.22.l2.ldscore.gz | output files from make_annot_from_geneset_all_chr.py follow the pattern <PREFIX>.<ANNOTATIONNAME>.<CHR>.l2.ldscore.gz


dict_genomic_annot = {}
ldscore_check_dict = collections.defaultdict(list) # checks that ldscores have been computed for all chromosomes for all cts
for file_ldscore in file_ldscores:
	m = re.search(r"^(.*)\.(\d{1,2})\.l2.ldscore.gz$", os.path.basename(file_ldscore)) # file_ldscore e.g. "/scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.yellow.9.l2.ldscore.gz"
	prefix_genomic_annot = m.groups()[0] # groups()[0]=prefix+annotation_name ; groups()[1]=chromosome
	chromosome = m.groups()[1]
	dict_genomic_annot[prefix_genomic_annot] = "{}/{}.".format(os.path.dirname(file_ldscore), prefix_genomic_annot) # *OBS*: dot is important
		# ^ prefix_genomic_annot = nn_lira_sema.yellow
		# ^ dict_genomic_annot[prefix_genomic_annot] =/scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.yellow.
	ldscore_check_dict[prefix_genomic_annot].append(int(chromosome))

### make check that annotation has ldscores for all chromosomes
chromosomes = set(range(1,23))
for prefix_genomic_annot in ldscore_check_dict.keys():
	ldscore_chromosomes = set(ldscore_check_dict[prefix_genomic_annot])
	if not ldscore_chromosomes == chromosomes: # 1...22
		chromosomes_missing = chromosomes.difference(ldscore_chromosomes)
		print("*WARNING*: prefix_genomic_annot={} did not have ldscore files for all chromosomes: {}. It will be dropped".format(prefix_genomic_annot, ",".join(map(str, chromosomes_missing))))
		dict_genomic_annot.pop(prefix_genomic_annot, None) # drop key from dict. REF: https://stackoverflow.com/a/11277439/6639640

### filter out 'all_genes' annotation
ANNOT_NAME_ALL_GENES="all_genes_in_dataset"
for prefix_genomic_annot in dict_genomic_annot.keys():
	# m = re.search(r".*%s.*" % ANNOT_NAME_ALL_GENES, os.path.basename(prefix_genomic_annot)) # REF 'using a variable inside a regex' https://stackoverflow.com/a/6931048/6639640
	if ANNOT_NAME_ALL_GENES in os.path.basename(prefix_genomic_annot): # check if string 'contains' substring
		print("*OBS*: prefix_genomic_annot={} matched ANNOT_NAME_ALL_GENES={}: It will be dropped to avoid colinearity in LDSC regression.".format(prefix_genomic_annot, ANNOT_NAME_ALL_GENES))
		dict_genomic_annot.pop(prefix_genomic_annot, None) # drop key from dict. REF: https://stackoverflow.com/a/11277439/6639640

with open(args.cts_outfile, "w") as fh_out:
	print("Writing output cts file: {}".format(args.cts_outfile))
	for prefix_genomic_annot in sorted(dict_genomic_annot):
		fh_out.write("{}\t{}\n".format(prefix_genomic_annot, dict_genomic_annot[prefix_genomic_annot]))


###################################### CALL ######################################





