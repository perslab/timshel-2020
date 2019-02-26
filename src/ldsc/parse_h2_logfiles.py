import os
import sys

import numpy as np
import pandas as pd

import glob
import re


def parse_ldsc_h2_logfile(filepath):
	name_gwas = re.search(r"(.*).log", os.path.basename(filepath)).group(1)
	print("Parsing gwas = {}".format(name_gwas))
	flag_extracted_data = False
	with open(filepath, "r") as fh_in:
		for line in fh_in:
			if line.startswith("Total Observed scale h2:"): # Total Observed scale h2: 0.1405 (0.0047)
				m = re.search(r"Total Observed scale h2: (.*) \((.*)\)", line)
				h2=float(m.group(1))
				h2_se=float(m.group(2))
				flag_extracted_data = True
			elif line.startswith("Total Liability scale h2:"): # Total Liability scale h2: 0.2164 (0.0081)
				m = re.search(r"Total Liability scale h2: (.*) \((.*)\)", line)
				h2=float(m.group(1))
				h2_se=float(m.group(2))
				flag_extracted_data = True
	if not flag_extracted_data:
		raise Exception("Did not find data in .log file: {}".format(filepath))
	return [name_gwas, h2, h2_se]
			

### Get files
dir_h2 = "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc_h2_baseline"
list_log_filepaths = glob.glob("{}/*.log".format(dir_h2))

### parse 
list_of_lists_parsed_data = []
for f in list_log_filepaths:
	list_of_lists_parsed_data.append(parse_ldsc_h2_logfile(f))
# list_of_lists_parsed_data[:5]


### Creating Pandas DataFrame from lists of lists.
df_h2 = pd.DataFrame(list_of_lists_parsed_data, columns=["gwas", "h2", "h2_se"])
df_h2.sort_values(by="gwas", inplace=True)
df_h2["h2_zscore"] = df_h2["h2"]/df_h2["h2_se"]
# df_h2.head()


### export
df_h2.to_csv("/projects/timshel/sc-genetics/sc-genetics/results/h2_observed_scale.multi_gwas.csv", index=False)


print("Script is done!")







