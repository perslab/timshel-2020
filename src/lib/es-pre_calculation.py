
# python3
import os
import sys
import pandas as pd
import numpy as np
from scipy import stats # f_oneway


def test_function(input_str):
	print(input_str)

def ctc_log_normalize(df):
	"""
	Perform common transcript count (ctc) normalization and log-transformation on input data
	
	Args:
		df:           DataFrame.
	"""
	print("Performning common transcript count (ctc) normalization and log-transformation on input data")
	### OBS INEFFICIENT: this creates a COPY of the data frame.
	df_ctc_log = np.log(1+df/df.sum(axis=0)*1e4) # column sum. Seurat default scale.factor is '10000'
	return(df_ctc_log)



def calculate_anova_sporadically_expressed_genes(df, annotations, out_prefix):
	"""
	See calculate_per_anno_summary_stats() for documentation
	"""
	print("Splitting data frame into annotation groups")
	annotations_unique = np.unique(annotations)
	list_arrays = []
	for counter, annotation in enumerate(annotations_unique, start=1):
		print("Splitting annotation #{}/#{} into group".format(counter, len(annotations_unique)))
		list_arrays.append(df.iloc[:, np.isin(annotations, [annotation])])
	### OLD - works
	# for lvl in df.columns.get_level_values(level="ClusterName").unique():
	# 	print(lvl)
	# 	list_arrays.append(df.xs(lvl, level="ClusterName", axis=1)) # # REF xs() : https://pandas.pydata.org/pandas-docs/stable/advanced.html#cross-section
	

	print("Running ANOVA")
	df_anova = pd.DataFrame(index=df.index)
	for idx_gene in range(len(df.index)): # loop over all genes
		if idx_gene % 100 == 0:
			print("gene {} out of {}".format(idx_gene, len(df.index)))
		gene_id = df.index[idx_gene]
		f = stats.f_oneway(*[df_tmp.values[idx_gene,:] for df_tmp in list_arrays]) # https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.stats.f_oneway.html
		df_anova.loc[gene_id, "pvalue"] = f.pvalue
		df_anova.loc[gene_id, "statistic"] = f.statistic
	df_anova.sort_values("pvalue",ascending=False, inplace=True)

	df_anova.to_csv("{}.pre_calc.sporadically_expressed_genes.anova.csv.gz".format(out_prefix), compression="gzip")

	n_genes_sporadically_expressed = np.sum(df_anova['pvalue'] > 0.00001) # 0.00001 is Skene cut-off
	print("Number of genes sporadically expressed (pvalue > 0.00001, Skene cut-off): {}".format(n_genes_sporadically_expressed))
	return(df_anova)



def calculate_per_anno_summary_stats(df, annotations, out_prefix, permute_annotations, seed=1):
	"""
	Pre-calculation of mean, variance and fraction expressed per annotation.
	
	Args:
		df:          		 DataFrame. genes x cells. Must have row index with gene IDs and column index. Any column indexes in the data frame will not be used (annotations argument defines the grouping of the data frame), so the DataFrame can have any kind of column index.
		annotations 		 1d array-like. Defines the annotations in df. Must have same length as number of columns in df.
		permute_annotations	 boolean. If true, annotations will be permuted for 'null generation'.
		seed: 	     		 int. Sets seed for annotation permutation. It is necessary to change the seed, if you want to generate different null backgrounds. Argument used if permute_annotations=T.
	Returns:
		DataFrame of mean, variance and fraction expressed. Each DataFrame has the same index (genes) and columns (unique annotations).
		DataFrame of number of cells per annotation. Single column ("n") and index is unique annotations. 
		Columns are alpha-numeric SORTED by their annotation name.
	Notes:
		Pandas works fine with non-unique column index. See some of the behavior here: https://github.com/pandas-dev/pandas/pull/3683
	"""
	### RUNTIME 'smart way' vectorized: <2 min
	### RUNTIME 'smart way' non-vectorized: 3h 8m 0s (~0.7 min per cell-type)
	### RUNTIME {self + others} 'dum way': ~5-7 min per cell-type (if doing mean and frac) --> ESTIMATED TOTAL RUNTIME = ~20 hours

	
	### Reorder data frame by sorted annotations. This should give a faster pandas data retrieval
	### Option 1: will duplicate object in memory, so slow and inefficient
	# df.columns = annotations # set column index. Duplicates are allowed.
	# df.sort_index(axis=1, inplace=True) # sort column index inplace
	### Option 2: alo slow and inefficien
	# np.argsort(): returns the indices that would sort an array.
	# df = df.iloc[:, np.argsort(annotations)] # COPIES the the data frame. Inefficient
	# REF: https://stackoverflow.com/a/39237712/6639640
	
	if permute_annotations:
		print("Doing null computation. Permuting labels with seed({}).".format(seed)
		np.random.seed(seed)
		annotations = np.random.permutation(annotations) # permute labels

	annotations_unique = np.unique(annotations) # returns SORTED unique elements of an array. We sort to make sure the output column is always in the same order
	
	df_frac = pd.DataFrame(index=df.index)
	df_mu = pd.DataFrame(index=df.index)
	df_var = pd.DataFrame(index=df.index)
	df_n = pd.DataFrame(index=annotations_unique) # obs different index

	df_frac.index.name="gene"
	df_mu.index.name="gene"
	df_var.index.name="gene"
	df_n.index.name="annotation"

	for counter, annotation in enumerate(annotations_unique, start=1):
		print("Running: #{}/#{} | {}".format(counter, len(annotations_unique), annotation))
		df_tmp_cells_in_annotation = df.iloc[:, np.isin(annotations, [annotation])] # boolean indexing. extract data once, then do computations. Data extraction is the slowest part of this function
		n_cells_in_annotation = df_tmp_cells_in_annotation.shape[1]
		df_n.loc[annotation, "n"] = n_cells_in_annotation
		df_frac.loc[:, annotation] = np.count_nonzero(df_tmp_cells_in_annotation, axis=1)/float(n_cells_in_annotation) # axis=1 : count non-zeros *along* columns | returns number of non-zeroes for each row.
		df_mu.loc[:, annotation] = df_tmp_cells_in_annotation.mean(axis='columns') # axis='columns': apply function to each row.
		df_var.loc[:, annotation] = df_tmp_cells_in_annotation.var(axis='columns') # Normalized by N-1 by default

	print("Writing files...")
	df_frac.to_csv("{}.pre_calc.frac_expr.{}csv.gz".format(out_prefix, "null." if permute_annotations else ""), compression="gzip")
	df_mu.to_csv("{}.pre_calc.mean.{}csv.gz".format(out_prefix, "null." if permute_annotations else ""), compression="gzip")
	df_var.to_csv("{}.pre_calc.var.{}csv.gz".format(out_prefix, "null." if permute_annotations else ""), compression="gzip")
	df_n.to_csv("{}.pre_calc.ncells.{}csv.gz".format(out_prefix, "null." if permute_annotations else ""), compression="gzip")

	return(df_frac, df_mu, df_var, df_n)

