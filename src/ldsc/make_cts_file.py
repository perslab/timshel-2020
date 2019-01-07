
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
	Will make cts files for ldscores files matching <prefix_ldscore_files>*.l2.ldscore.gz UNLESS --annotation_filter given.
	Notice that if prefix_ldscore_files is a directory (with a trailing '/'), then all ldscore file in the directory will be used.
	Example with prefix as input: /scratch/sc-ldsc/ldsc_analysis_dir/files_prefix
	Example with dir as input: /scratch/sc-ldsc/ldsc_analysis_dir/ [OBS the trailing slash is important!]""")
parser.add_argument('--cts_outfile', type=str, help="Filename for output CTS file.", required=True)
parser.add_argument('--annotation_filter', type=str, 
	help="""If specified, the cts_outfile will only contain annotations matching annotation names given in the annotation_filter file. 
	Argument value should be a filename to a file containing 'annotation filters'. 
	File should contain one annotation name per line for the annotations. 
	Annotations are filtered by including only those .l2.ldscore.gz files that contains the annotation filter substring.""")
# REF required=True: https://docs.python.org/3/library/argparse.html#required

args = parser.parse_args()


file_ldscores = glob.glob("{}*.l2.ldscore.gz".format(args.prefix_ldscore_files)) # e.g. <SOMEPATH>/nn_lira_sema.yellow.22.l2.ldscore.gz | output files from make_annot_from_geneset_all_chr.py follow the pattern <PREFIX>.<ANNOTATIONNAME>.<CHR>.l2.ldscore.gz


dict_genomic_annot = {}
ldscore_check_dict = collections.defaultdict(list) # checks that ldscores have been computed for all chromosomes for all cts
for file_ldscore in file_ldscores:
	m = re.search(r"^(.*)\.(\d{1,2})\.l2.ldscore.gz$", os.path.basename(file_ldscore)) # file_ldscore e.g. "/scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.yellow.9.l2.ldscore.gz"
	basename_prefix_genomic_annot_and_annot_name = m.groups()[0] # groups()[0]=prefix+annotation_name ; groups()[1]=chromosome
	# TODO POTENTIAL: to extract only annot_name, split basename_prefix_genomic_annot_and_annot_name on double underscore (if name contains it).
	chromosome = m.groups()[1]
	dict_genomic_annot[basename_prefix_genomic_annot_and_annot_name] = "{}/{}.".format(os.path.dirname(file_ldscore), basename_prefix_genomic_annot_and_annot_name) # *OBS*: dot is important
		# ^ basename_prefix_genomic_annot_and_annot_name = nn_lira_sema.yellow
		# ^ dict_genomic_annot[basename_prefix_genomic_annot_and_annot_name] =/scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.yellow.
	ldscore_check_dict[basename_prefix_genomic_annot_and_annot_name].append(int(chromosome))
print("Found n={} annotations in --prefix_ldscore_files".format(len(dict_genomic_annot)))

### make check that annotation has ldscores for all chromosomes
chromosomes = set(range(1,23))
for basename_prefix_genomic_annot_and_annot_name in list(ldscore_check_dict.keys()): # list() needed for py3 compatibility
	ldscore_chromosomes = set(ldscore_check_dict[basename_prefix_genomic_annot_and_annot_name])
	if not ldscore_chromosomes == chromosomes: # 1...22
		chromosomes_missing = chromosomes.difference(ldscore_chromosomes)
		print("*WARNING*: basename_prefix_genomic_annot_and_annot_name={} did not have ldscore files for all chromosomes: {}. It will be dropped".format(basename_prefix_genomic_annot_and_annot_name, ",".join(map(str, chromosomes_missing))))
		dict_genomic_annot.pop(basename_prefix_genomic_annot_and_annot_name, None) # drop key from dict while iterating over it. REF: https://stackoverflow.com/questions/5384914/how-to-delete-items-from-a-dictionary-while-iterating-over-it and https://stackoverflow.com/a/11277439/6639640

# ### Filter out 'all_genes' annotation [***LEGACY CODE*** - CAN BE REMOVED LATER. all_genes annotations are new kept in seperate directory]
# ANNOT_NAME_ALL_GENES="all_genes_in_dataset"
# for basename_prefix_genomic_annot_and_annot_name in list(dict_genomic_annot.keys()): # list() needed for py3 compatibility
# 	# m = re.search(r".*%s.*" % ANNOT_NAME_ALL_GENES, os.path.basename(basename_prefix_genomic_annot_and_annot_name)) # REF 'using a variable inside a regex' https://stackoverflow.com/a/6931048/6639640
# 	if ANNOT_NAME_ALL_GENES in os.path.basename(basename_prefix_genomic_annot_and_annot_name): # check if string 'contains' substring
# 		print("*OBS*: basename_prefix_genomic_annot_and_annot_name={} matched ANNOT_NAME_ALL_GENES={}: It will be dropped to avoid colinearity in LDSC regression.".format(basename_prefix_genomic_annot_and_annot_name, ANNOT_NAME_ALL_GENES))
# 		dict_genomic_annot.pop(basename_prefix_genomic_annot_and_annot_name, None) # drop key from dict while iterating over it. REF: https://stackoverflow.com/a/11277439/6639640


### Apply annotation filter
if args.annotation_filter:
	print("Applying annotation file using file: {}".format(args.annotation_filter))
	lines_annotation_filter = open(args.annotation_filter, 'r').read().splitlines() # read file and split on line breaks. splitlines() handles universal newline. REF: https://stackoverflow.com/a/12330535/6639640
	print("Read n={} annotation filters.".format(len(lines_annotation_filter)))
	for basename_prefix_genomic_annot_and_annot_name in list(dict_genomic_annot.keys()): # list() needed for py3 compatibility
		if not any(filter_name in os.path.basename(basename_prefix_genomic_annot_and_annot_name) for filter_name in lines_annotation_filter): # check if any of the annotation_filters names are contained in basename_prefix_genomic_annot_and_annot_name as a substring. REF: https://stackoverflow.com/a/3389611/6639640
			dict_genomic_annot.pop(basename_prefix_genomic_annot_and_annot_name, None) # drop key from dict while iterating over it. REF: https://stackoverflow.com/questions/5384914/how-to-delete-items-from-a-dictionary-while-iterating-over-it and https://stackoverflow.com/a/11277439/6639640
print("After filtering, n={} annotations are retained".format(len(dict_genomic_annot)))

if len(dict_genomic_annot) == 0:
	raise Exception("dict_genomic_annot is empty and no annotations can be written to cts outfile. Check that your annotation filters are correct.")

with open(args.cts_outfile, "w") as fh_out:
	print("Writing output cts file: {}".format(args.cts_outfile))
	for basename_prefix_genomic_annot_and_annot_name in sorted(dict_genomic_annot):
		fh_out.write("{}\t{}\n".format(basename_prefix_genomic_annot_and_annot_name, dict_genomic_annot[basename_prefix_genomic_annot_and_annot_name]))


###################################### CALL ######################################





