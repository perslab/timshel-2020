#!/usr/bin/env python2.7
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
# import gzip


### PT adds
import multiprocessing
import os


###################################### USAGE ######################################
# Must run on python2.7, because of issue with pybedtools in python3.

### test
# python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /raid5/projects/timshel/sc-genetics/ldsc/ldsc/test_file_multi_gene_set_wgcna200.csv \
# --file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir  ./tmp_tmp_test_make_annot \
# --out_prefix test_xxx1 \
# --flag_wgcna \
# --flag_mouse


### lira sema (scratch)
# python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /raid5/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv \
# --file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir  /scratch/sc-ldsc/nn_lira_sema \
# --out_prefix nn_lira_sema \
# --flag_wgcna \
# --flag_mouse



###################################### FILE SNIPPETS ######################################


# @@@@@@@@@@@@@@@@@@ file_multi_gene_set @@@@@@@@@@@@@@@@@@

### WGCNA input (flag_wgcna = True). 
### We use the 'module' column as annotation name, since they should be unique
# cell_cluster,module,ensembl,hgcn,pkME
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000017639,Rab11fip4,0.888954070670456
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000004233,Wars2,0.860974901276187
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000031487,Brf2,0.855919687669775
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000032997,Chpf,0.824264829666081

### non WGCNA input.
### No header. Delim=csv. 
### Col1=annotation_name, Col2=EnsemblID (mouse or human)
# brain_cortex,ENSMUSG00000017639
# brain_cortex,ENSMUSG00000004233
# brain_cortex,ENSMUSG00000031487
# brain_cortex,ENSMUSG00000032997
# ...
# brain_hypothlamus,ENSMUSG00000032997

def map_ensembl_genes_mouse_to_human(df):
    """
    DESCRIPTION
        function to map mouse ensembl genes to human ortholog genes.
    INPUT
       df: a dataframe with two columns: "annotation" and "gene". "gene" column should contain Ensembl mouse gene IDs to be mapped.
    OUTPUT
       df: input dataframe with mapped genes. Mouse genes that could not be mapped are removed.
       file_mapping_summary: a summary file with mapping stats

    REMARKS
        We assume the file_mapping contains only 1-1 mapping (which is true for gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz).
        Otherwise the .map() function might fail
    """
    file_out_mapping_stats = "{}.make_annotation_mapping_stats.txt".format(out_prefix)
    
    file_mapping = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
    ### SNIPPET
    # ensembl_gene_id chromosome_name start_position  end_position    mmusculus_homolog_ensembl_gene  mmusculus_homolog_orthology_confidence
    # ENSG00000138593 15      49280673        49338760        ENSMUSG00000035093      1
    # ENSG00000166351 21      14982498        15013906        ENSMUSG00000095294      0
    # ENSG00000168675 18      13217497        13652754        ENSMUSG00000024544      1
    
    df_mapping = pd.read_csv(file_mapping, delim_whitespace=True, usecols=["ensembl_gene_id", "mmusculus_homolog_ensembl_gene"], index_col=1) # index is MOUSE genes
    
    ### map ortholog genes
    df["gene_human"] = df["gene"].map(df_mapping["ensembl_gene_id"]) # df["gene"] is mouse genes. df_mapping["ensembl_gene_id"] is human.
    # ^ .map() returns NaN values for genes not mapped. 

    ### make summary of mapping
    df_summary = df.groupby("annotation")["gene_human"].agg({'n_genes_input': lambda x: len(x),
                                                         'n_genes_output': lambda x: len(x)-sum(pd.isnull(x)), 
                                                        'n_genes_not_mapped' : lambda x: sum(pd.isnull(x)),
                                                        'pct_genes_not_mapped': lambda x: "{:.2f}".format(sum(pd.isnull(x))/float(len(x))*100)})
    # ^ we use pd.isnull() instead of np.isnan() because of the issue described here: https://stackoverflow.com/a/36001191/6639640
    df_summary.sort_values(by=['n_genes_output'], inplace=True) 
    df_summary.to_csv(file_out_mapping_stats, sep="\t")
    
    ### final processing
    df.dropna(axis=0, inplace=True) # remove non-mapped genes
    df['gene'] = df['gene_human'] # replace mouse genes with human genes
    df.drop(['gene_human'], axis=1, inplace=True) # remove human gene column
    
    return df


def read_multi_gene_set_file(file_multi_gene_set, flag_wgcna, flag_mouse):
    if flag_wgcna:
        df_multi_gene_set = pd.read_csv(file_multi_gene_set, usecols=["module", "ensembl"]) 
        # ^ 'module' = first column; will later be renamed to 'annotation'
        # ^ 'ensembl' = second column; will later be renamed to 'gene'
    else:
        df_multi_gene_set = pd.read_csv(file_multi_gene_set, header=None)
    df_multi_gene_set.columns = ["annotation", "gene"] # setting or renaming column names
    if flag_mouse:
        df_multi_gene_set = map_ensembl_genes_mouse_to_human(df_multi_gene_set)
    return df_multi_gene_set


def multi_gene_sets_to_dict_of_beds(df_multi_gene_set, df_gene_coord, windowsize):
    """ 
    INPUT
        df_multi_gene_set: two columns "annotation" and "gene". Gene is human Ensembl gene names.
    OUTPUT
        dict_of_beds: returns a dict of beds. Keys are annotation names from df_multi_gene_set.
    """
    print('making gene set bed files')
    dict_of_beds = {}
    for name_annotation, df_group in df_multi_gene_set.groupby("annotation"):
        df = pd.merge(df_gene_coord, df_group, left_on="GENE", right_on="gene", how = "inner") 
        df['START'] = np.maximum(0, df['START'] - windowsize)
        df['END'] = df['END'] + windowsize
        iter_df = [['chr'+(str(x1).lstrip('chr')), x2, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])] # notice that only 3 columns from df is used.
        bed_for_annot = BedTool(iter_df).sort().merge()
        dict_of_beds[name_annotation] = bed_for_annot
    return dict_of_beds

def make_annot_file_one_chromosome(chromosome):
    """ 
    Input
        chromosome: integer (1..22)
    
    *OBS* this function RELIES on MANY GLOBAL scope VARIABLES
    """
    # TODO: parse variables to function
    print('making annot files for chromosome {}'.format(chromosome))
    bimfile = "{}.{}.bim".format(args.bimfile_basename, chromosome) # e.g. <LONGPATH/1000G.EUR.QC>.15.bim for chromosome 15
    df_bim = pd.read_csv(bimfile, delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [['chr'+str(x1), x2, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)
    counter = 1
    for name_annotation in dict_of_beds:
        print("CHR={} | annotation={}, #{}/#{}".format(chromosome, name_annotation, counter, len(dict_of_beds)))
        file_out_annot = "{}/{}.{}.{}.annot.gz".format(out_dir, out_prefix, name_annotation, chromosome) # set output filename. ${prefix}.${chr}.annot.gz
        bed_for_annot = dict_of_beds[name_annotation] # get bed
        
        annotbed = bimbed.intersect(bed_for_annot)
        bp = [x.start for x in annotbed]
        df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
        df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
        df_annot.fillna(0, inplace=True)
        df_annot = df_annot[['ANNOT']].astype(int)
        df_annot.to_csv(file_out_annot, sep="\t", index=False, compression="gzip")
        counter += 1
    return None



###################################### MAIN ######################################

parser = argparse.ArgumentParser()
parser.add_argument('--file_multi_gene_set', type=str, help='a file of gene names, one line per gene.')
parser.add_argument('--file_gene_coord', type=str, help='a file with columns GENE, CHR, START, and END, where START and END are base pair coordinates of TSS and TES. This file can contain more genes than are in the gene set. We provide ENSG_coord.txt as a default.')
parser.add_argument('--windowsize', type=int, help='how many base pairs to add around the transcribed region to make the annotation? Finucane uses 100 kb.')
parser.add_argument('--bimfile_basename', type=str, help='plink bim BASENAME for the dataset you will use to compute LD scores. If argument is "1000G.EUR.QC", then the files "1000G.EUR.QC.1.bim", "1000G.EUR.QC.2.bim", ..., "1000G.EUR.QC.22.bim" will be loaded')
parser.add_argument('--out_dir', type=str, help='output directory to write annot files. Relative or absolute. Dir be created if it does not exist. ')
parser.add_argument('--out_prefix', type=str, help='Prefix for output files. Outputfiles will be <out_dir>/<out_prefix>.<name_annotation>.<chromosome>.annot.gz')
parser.add_argument('--flag_wgcna', action='store_true', help='set flag if file_multi_gene_set input is from WGCNA pipeline')
parser.add_argument('--flag_mouse', action='store_true', help='set flag if ile_multi_gene_set input has mouse genes (instead of human)')
# parser.add_argument('--n_proc', type=int, help='Number of processes. Default 22 (the number of chromosomes)')

#parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')


args = parser.parse_args()

### TODO
# [ ] check that ALL bim files exist. (so pool of workers does not fail)

out_prefix = args.out_prefix
out_dir = args.out_dir

### Make out_dir
if not os.path.exists(out_dir):
    print("Making output dir {}".format(out_dir))
    os.makedirs(out_dir)


### read coord file
df_gene_coord = pd.read_csv(args.file_gene_coord, delim_whitespace = True)
### read gene list file (a list of gene list files)
df_multi_gene_set = read_multi_gene_set_file(args.file_multi_gene_set, flag_wgcna=args.flag_wgcna, flag_mouse=args.flag_mouse)
### make beds
dict_of_beds = multi_gene_sets_to_dict_of_beds(df_multi_gene_set, df_gene_coord, args.windowsize)

print("Starting pool...")
pool = multiprocessing.Pool(processes=22)
pool.map(make_annot_file_one_chromosome, range(1,23)) # 1...22

print("Script is done!")



