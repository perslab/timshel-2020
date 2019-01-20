
import glob
import re
import os



models = ["baseline.all_genes",
"baseline.no_all_genes",
"no_baseline.all_genes",
"no_baseline.no_all_genes"]

header = ["Model", # PT
"Annotation", # PT
"Category",
"Prop._SNPs",
"Prop._h2",
"Prop._h2_std_error",
"Enrichment",
"Enrichment_std_error",
"Enrichment_p",
"Coefficient",
"Coefficient_std_error",
"Coefficient_z-score"]

all_model_results = []
for model in models:
	files = glob.glob("/raid5/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/experiment-h2.{model}/*.results".format(model=model))
	model_results = []
	for file_result in files:
		annotation_name = re.search(r"experiment.gtex_tissues.h2.(.*).BMI_Yengo2018.results$", os.path.basename(file_result)).groups()[0] # experiment.gtex_tissues.h2.Artery_Coronary.BMI_Yengo2018.results
		# print(annotation_name)
		with open(file_result, "r") as f:
			result = f.readlines()[1].strip() # 2nd line in file contain our result
			result_line = "{}\t{}\t{}".format(model, annotation_name, result)
			model_results.append(result_line)
	model_outfile = "experiment.gtex_tissues.h2.{model}.BMI_Yengo2018.txt".format(model=model)
	with open(model_outfile, "w") as fout:
		fout.write("\t".join(header)+"\n")
		for l in model_results:
			fout.write(l+"\n")
	all_model_results.extend(model_results)


# all model results
all_model_outfile = "experiment.gtex_tissues.h2.{model}.BMI_Yengo2018.txt".format(model="ALL_MODELS")
with open(all_model_outfile, "w") as fout:
	fout.write("\t".join(header)+"\n")
	for l in all_model_results:
		fout.write(l+"\n")




