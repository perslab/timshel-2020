from __future__ import print_function

import multiprocessing
import os
import sys
import subprocess

import itertools
import functools

import pandas as pd
import numpy as np
import argparse

from pybedtools import BedTool 
# ^ some/all of pybedtools requires 'bedtools' to be available on your PATH /tools/bedtools/2.27.1/bin/bedtools
# BedTool(..).sort(): "sortBed" does not appear to be installed or on the path, so this method is disabled. Please install a more recent version of BEDTools and re-import to use this method.
# Install via conda install --channel conda-forge --channel bioconda pybedtools bedtools htslib
# [did not work for me, so I used conda install ... instead] Use the below lines to add a system installation of bedtools and tabix when running script within anaconda:
# os.environ["PATH"] += os.pathsep + "/tools/bedtools/2.27.1/bin/"
# os.environ["PATH"] += os.pathsep + "/tools/htslib/1.6/bin/"

import psutil

import pdb

# import gzip




###################################### USAGE ######################################

### test
# python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /raid5/projects/timshel/sc-genetics/ldsc/ldsc/test_file_multi_gene_set_wgcna200.csv \
# --file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir ./tmp_tmp_test_make_annot \
# --out_prefix test_xxx1 \
# --flag_wgcna \
# --flag_mouse


### lira sema (scratch)
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /raid5/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv \
# --file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/nn_lira_sema \
# --out_prefix nn_lira_sema \
# --flag_wgcna \
# --flag_mouse


### MACA
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/jonatan/tmp-maca/tables/maca_tissue_cell_type_kME_cell_cluster_module_genes.csv \
# --file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/maca \
# --out_prefix maca_tissue_cell_type \
# --flag_wgcna \
# --flag_mouse

### Mousebrain Neurons_ClusterName (n=3833)
### 11 processes takes ~500 GB
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/jonatan/tmp-mousebrain/tables/mousebrain_Neurons_ClusterName_2_cell_cluster_module_genes.csv \
# --file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/mousebrain \
# --out_prefix Neurons_ClusterName \
# --flag_wgcna \
# --flag_mouse \
# --n_parallel_jobs 3


### Mette thesis hypothalamus
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /raid5/projects/timshel/sc-genetics/sc-genetics/data/gene_lists/mludwig_thesis_hypothalamus_wgcna_modules.csv \
# --file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/hypothalamus_mette_thesis \
# --out_prefix hypothalamus_mette_thesis \
# --flag_wgcna \
# --flag_mouse


### Mousebrain Ependymal_ClusterName (n=XXX)
# time python2 make_annot_from_geneset_all_chr.py \
# --file_multi_gene_set /projects/jonatan/tmp-mousebrain/tables/mousebrain_Ependymal_ClusterName_2_cell_cluster_module_genes.csv \
# --file_gene_coord /raid5/projects/timshel/sc-genetics/ldsc/data/gene_coords/gene_annotation.hsapiens_all_genes.GRCh37.ens_v91.LDSC_fmt.txt \
# --windowsize 100000 \
# --bimfile_basename /raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC \
# --out_dir /scratch/sc-ldsc/mousebrain \
# --out_prefix Ependymal_ClusterName \
# --flag_wgcna \
# --flag_mouse \
# --n_parallel_jobs 3



###################################### DOCUMENTATION ######################################

### General remarks:
# - any whitespace in annotation_name column in file_multi_gene_set will be converted to underscore ('_'). 
#   This is because LDSC .annot files are read as *whitespace delimted* by the ldsc.py program, so annotation_name with whitespace in the name will make the .l2.ldscore.gz header wrong.

### For continuous annotations:
# - When a variant is spanned by multiple genes with the XXX kb window, we assign the maximum annotation_value.

### For binary annotations:
# - The SNPs within the genomic regions spanned by the genes within a given annotation gets the annotation value 1. All other SNPs get the annotation value 0



###################################### FILE SNIPPETS ######################################


# @@@@@@@@@@@@@@@@@@ file_multi_gene_set @@@@@@@@@@@@@@@@@@

### WGCNA input (flag_wgcna = True). 
### We use the 'module' column as annotation name, since they should be unique
### See read_multi_gene_set_file() for what columns are used in the file.
# cell_cluster,module,ensembl,hgcn,pkME
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000017639,Rab11fip4,0.888954070670456
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000004233,Wars2,0.860974901276187
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000031487,Brf2,0.855919687669775
# Aorta_endothelial cell,antiquewhite3,ENSMUSG00000032997,Chpf,0.824264829666081

### non WGCNA input.
### NO HEADER. Delim=csv. 
### Col1=annotation_name, Col2=EnsemblID (mouse or human), Col3
### Col3 is not required (and not used) if --flag_encode_as_binary_annotation is set
# brain_cortex,ENSMUSG00000017639,0.34
# brain_cortex,ENSMUSG00000004233,0.01
# brain_cortex,ENSMUSG00000031487,0.22
# brain_cortex,ENSMUSG00000032997,0.98
# ...
# brain_hypothlamus,ENSMUSG00000032997, 0.87


def map_ensembl_genes_mouse_to_human(df, args):
    """
    DESCRIPTION
        function to map mouse ensembl genes to human ortholog genes.
    INPUT
       df: a dataframe with two columns: "annotation" and "gene_input". "gene_input" column should contain Ensembl mouse gene_input IDs to be mapped.
    OUTPUT
       df: input dataframe with mapped genes. Mouse genes that could not be mapped are removed.
       file_mapping_summary: a summary file with mapping stats

    REMARKS
        We assume the file_mapping contains only 1-1 mapping (which is true for gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz).
        Otherwise the .map() function might fail
    """
    print("Mapping from mouse to human genes")
    file_out_mapping_stats = "{}/log.{}.make_annotation_mapping_stats.txt".format(args.out_dir, args.out_prefix)
    
    file_mapping = "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
    ### SNIPPET
    # ensembl_gene_id chromosome_name start_position  end_position    mmusculus_homolog_ensembl_gene  mmusculus_homolog_orthology_confidence
    # ENSG00000138593 15      49280673        49338760        ENSMUSG00000035093      1
    # ENSG00000166351 21      14982498        15013906        ENSMUSG00000095294      0
    # ENSG00000168675 18      13217497        13652754        ENSMUSG00000024544      1
    
    df_mapping = pd.read_csv(file_mapping, delim_whitespace=True, usecols=["ensembl_gene_id", "mmusculus_homolog_ensembl_gene"], index_col=1) # index is MOUSE genes
    
    ### map ortholog genes
    df["gene"] = df["gene_input"].map(df_mapping["ensembl_gene_id"]) # df["gene_input"] is mouse genes. df_mapping["ensembl_gene_id"] is human.
    # ^ .map() returns NaN values for genes not mapped. 

    ### make summary of mapping
    df_summary = df.groupby("annotation")["gene"].agg({'n_genes_input': lambda x: len(x),
                                                         'n_genes_output': lambda x: len(x)-sum(pd.isnull(x)), 
                                                        'n_genes_not_mapped' : lambda x: sum(pd.isnull(x)),
                                                        'pct_genes_not_mapped': lambda x: "{:.2f}".format(sum(pd.isnull(x))/float(len(x))*100)})
    # ^ we use pd.isnull() instead of np.isnan() because of the issue described here: https://stackoverflow.com/a/36001191/6639640
    df_summary.sort_values(by=['n_genes_output'], inplace=True) 
    df_summary.to_csv(file_out_mapping_stats, sep="\t")
    
    ### final processing
    df.dropna(axis=0, inplace=True) # remove non-mapped genes
    return df


def read_multi_gene_set_file(args):
    file_multi_gene_set = args.file_multi_gene_set
    if args.flag_wgcna:
        df_multi_gene_set = pd.read_csv(file_multi_gene_set) 
        df_multi_gene_set.rename(columns={"module":"annotation", "ensembl":"gene_input", "pkME":"annotation_value"}, inplace=True) # df.rename(columns={'oldName1': 'newName1', 'oldName2': 'newName2'})
        # ^ 'module' = first column; will later be renamed to 'annotation'
        # ^ 'ensembl' = second column; will later be renamed to 'gene'
    else:
        df_multi_gene_set = pd.read_csv(file_multi_gene_set, sep=None, header=None, engine='python') # sep=None: automatically detect the separator
        df_multi_gene_set.columns = df_multi_gene_set.columns.map(str) # because header=None, the .columns is of type integer (Int64Index(.., dtype='int64')). We need to map array to string before we can rename the columns below
        if args.flag_encode_as_binary_annotation:
            df_multi_gene_set.columns.values[[0,1]] = ["annotation", "gene_input"] # .values() is needed to avoid TypeError
            # ALTERNATIVE ---> https://stackoverflow.com/a/43759994/6639640: df.rename(columns={ df.columns[1]: "your value" })
        else:
            df_multi_gene_set.columns.values[[0,1,2]] = ["annotation", "gene_input", "annotation_value"]
    if args.flag_encode_as_binary_annotation:
        print("Converting to binary encoding")
        df_multi_gene_set["annotation_value"] = 1 # binary encoding. All genes get the value 1.
    else:
        if not np.issubdtype(df_multi_gene_set["annotation_value"].dtype, np.number): # REF: https://stackoverflow.com/a/38185759/6639640
            raise Exception("ERROR: your df_multi_gene_set contains non-numeric annotation values. Will not create annotation files.")
        if (df_multi_gene_set["annotation_value"] < 0).any():
            raise Exception("ERROR: your df_multi_gene_set contains negative annotation values. Will not create annotation files.")
    ### Make filepath and LDSC 'valid' annotation name (read: clean up the names!)
    df_multi_gene_set["annotation"] = df_multi_gene_set["annotation"].replace(r"\s+", "_",regex=True) 
    # ^ any whitespace in annotation_name column in file_multi_gene_set will be converted to underscore ('_').  
    # ^ This is because LDSC .annot files are read as *whitespace delimted* by the ldsc.py program, so annotation_name with whitespace in the name will make the .l2.ldscore.gz header wrong.
    df_multi_gene_set["annotation"] = df_multi_gene_set["annotation"].replace(r"/", "-",regex=True) 
    # ^ We need to avoid forward slash (/) in the filenames when files are split per annotation (/per_annot dir). 
    # ^ If forward slashes are not replaced, we would get an error when looking for or writing files such as "celltypes.campbell_lvl1.all.campbell_lvl1.a06.NG2/OPC.ges.21.l2.M" when the annotation name is "a06.NG2/OPC.ges"
    print("Read file_multi_gene_set. Header of the parsed/processed file:")
    print(df_multi_gene_set.head(10))
    print("Annotation value summary stats:")
    df_annot_value_sumstats = df_multi_gene_set.groupby("annotation")["annotation_value"].agg(["mean", "std", "max", "min", "count"])
    print(df_annot_value_sumstats)
    file_out_annot_value_sumstatsstats = "{}/log.{}.make_annotation_value_sumstats.txt".format(args.out_dir, args.out_prefix)
    df_annot_value_sumstats.to_csv(file_out_annot_value_sumstatsstats, sep="\t")


    print("========================== STATS file_multi_gene_set ====================")
    print("Number of gene sets: {}".format(df_multi_gene_set["annotation"].nunique()))
    print("=========================================================================")
    if args.flag_mouse:
        df_multi_gene_set = map_ensembl_genes_mouse_to_human(df_multi_gene_set, args) # adds "gene" column
    else:
        df_multi_gene_set["gene"] = df_multi_gene_set["gene_input"] # copy
    ### write_multi_geneset_file
    file_out_multi_geneset = "{}/log.{}.multi_geneset.txt".format(args.out_dir, args.out_prefix)
    df_multi_gene_set.to_csv(file_out_multi_geneset, sep="\t", index=False)
    return df_multi_gene_set


def multi_gene_sets_to_dict_of_beds(df_multi_gene_set, df_gene_coord, windowsize):
    """ 
    INPUT
        df_multi_gene_set: three columns "annotation", "gene" and "annotation_value". Gene is human Ensembl gene names.
    OUTPUT
        dict_of_beds: returns a dict of beds. Keys are annotation names from df_multi_gene_set.
    """
    print('making gene set bed files')
    #n_genes_not_in_gene_coord = np.sum(np.isin(df_multi_gene_set["gene"], df_gene_coord["GENE"], invert=True)) # numpy.isin(element, test_elements). Calculates element in test_elements, broadcasting over element only. Returns a boolean array of the same shape as element that is True where an element of element is in test_elements and False otherwise.
    #if n_genes_not_in_gene_coord > 0:
    #    print("*WARNING*: {} genes in the (mapped) input multi gene set is not found in the gene coordinate file. These genes will be discarded".format(n_genes_not_in_gene_coord))
    dict_of_beds = {}
    for name_annotation, df_group in df_multi_gene_set.groupby("annotation"):
        print("Merging input multi gene set with gene coordinates for annotation = {}".format(name_annotation))
        df = pd.merge(df_gene_coord, df_group, left_on="GENE", right_on="gene", how = "inner")
        df['START'] = np.maximum(0, df['START'] - windowsize)
        df['END'] = df['END'] + windowsize
        list_of_lists = [['chr'+(str(chrom).lstrip('chr')), str(start), str(end), str(name), str(score)] for (chrom,start,end,name,score) in np.array(df[['CHR', 'START', 'END', 'GENE', 'annotation_value']])]
        # consider using BedTool.from_dataframe(df[, outfile, sep, header, .])   Creates a BedTool from a pandas.DataFrame.
        # BedTool() can accept a list or tubple for creation: https://github.com/daler/pybedtools/blob/master/pybedtools/bedtool.py#L511
        bed_for_annot = BedTool(list_of_lists).sort().merge(c=[4,5], o=["distinct","max"]) 
            # ^ .merge(c=[5], o=["max"]): When a variant is spanned by multiple genes with the XXX kb window, we assign the maximum annotation_value (column 5, score field).
            # ^ .merge(c=[4], o=["distinct"]): bed_for_annot's name field (column 4) will be the distinct genes for in that feature. Fetures that was merged will contain multiple genes in the field.
            # ^ .merge(): Merge overlapping features together. https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.merge.html#pybedtools.bedtool.BedTool.merge
            # ^ .merge(): OBS: requires that you PRESORT your data by chromosome and then by start position
            # ^ .merge(): OBS: without any -c / -o flags, only the first 3 columns of the bed file will be kept (chrom, start, end)
            # ^ .merge(): both of these syntaxes should work: .merge(c="4,5", o="distinct,max") and .merge(c=[4,5], o=["distinct","max"]) 
        ### BED5 format. REF: https://bedtools.readthedocs.io/en/latest/content/general-usage.html
        # chrom - The name of the chromosome on which the genome feature exists.
        #     Any string can be used. For example, "chr1", "III", "myChrom", "contig1112.23".
        #     This column is required.
        # start - The zero-based starting position of the feature in the chromosome.
        #     The first base in a chromosome is numbered 0.
        #     The start position in each BED feature is therefore interpreted to be 1 greater than the start position listed in the feature. For example, start=9, end=20 is interpreted to span bases 10 through 20,inclusive.
        #     This column is required.
        # end - The one-based ending position of the feature in the chromosome.
        #     The end position in each BED feature is one-based. See example above.
        #     This column is required.
        # name - Defines the name of the BED feature.
        #     Any string can be used. For example, "LINE", "Exon3", "HWIEAS_0001:3:1:0:266#0/1", or "my_Feature".
        #     This column is optional.
        # score - The UCSC definition requires that a BED score range from 0 to 1000, inclusive. However, bedtools allows any string to be stored in this field in order to allow greater flexibility in annotation features. For example, strings allow scientific notation for p-values, mean enrichment values, etc. It should be noted that this flexibility could prevent such annotations from being correctly displayed on the UCSC browser.
        #     Any string can be used. For example, 7.31E-05 (p-value), 0.33456 (mean enrichment value), "up", "down", etc.
        #     This column is optional.
        dict_of_beds[name_annotation] = bed_for_annot
    return dict_of_beds


def get_annot_file_path(chromosome, args):
    file_out_annot_combined = "{}/{}.{}.{}.annot.gz".format(args.out_dir, args.out_prefix, "COMBINED_ANNOT", chromosome) # set output filename. ${prefix}.${chr}.annot.gz
    return file_out_annot_combined
    


def check_annot_file(chromosome, args):
    ### check for existing annot file (and it's integrity)
    file_out_annot_combined = get_annot_file_path(chromosome, args)
    if os.path.exists(file_out_annot_combined):
        print("CHR={} | file_out_annot_combined exists: {}. Will check for integrity...".format(chromosome, file_out_annot_combined))
        cmd_gzip_test = "gzip -t {}".format(file_out_annot_combined)
        with open(os.devnull, 'w') as fnull:
            p = subprocess.Popen(cmd_gzip_test, shell=True, stdout=fnull, stderr=subprocess.STDOUT)
            # gzip -t doesn't have any output, other than the return code, if it's a correct gzip compressed file. But if it is a corrupted file, it will output something:
            # CORRUPT FILE --> gzip -t Neurons_ClusterName.COMBINED_ANNOT.17.annot.gz --> "gzip: Neurons_ClusterName.COMBINED_ANNOT.17.annot.gz: unexpected end of file"
            p.wait()
        if p.returncode == 0: # gz file ok
            print("CHR={} | file_out_annot_combined integrity OK: {}. Will not make new file".format(chromosome, file_out_annot_combined))
            return None
        else: # file not ok.
            return chromosome
    else: # file does not exist
        return chromosome



def make_annot_file_per_chromosome(chromosome, dict_of_beds, args):
    """ 
    Input
        chromosome: integer (1..22)
    
    *OBS* this function RELIES on MANY GLOBAL scope VARIABLES
    """
    # TODO: parse variables to function
    
    ### make annot file
    print('making annot files for chromosome {}'.format(chromosome))
    bimfile = "{}.{}.bim".format(args.bimfile_basename, chromosome) # e.g. <LONGPATH/1000G.EUR.QC>.15.bim for chromosome 15
    df_bim = pd.read_csv(bimfile, delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    # (Pdb) df_bim.head()
    #    CHR          SNP        CM       BP
    # 0   21  rs146134162 -0.908263  9412099
    # 1   21  rs578050168 -0.908090  9412377
    # 2   21  rs527616997 -0.907297  9413645
    # 3   21  rs544748596 -0.906578  9414796
    # 4   21  rs528236937 -0.906500  9414921
    
    # iter_bim = [['chr'+str(x1), x2, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    # ^ Python3 (but not Python2.7) gives the following error when calling "bimbed = BedTool(iter_bim)" in make_annot_file_per_chromosome()
    # /tools/anaconda/3-4.4.0/lib/python3.6/site-packages/pybedtools/cbedtools.pyx in pybedtools.cbedtools.IntervalIterator.__next__()
    # /tools/anaconda/3-4.4.0/lib/python3.6/site-packages/pybedtools/cbedtools.pyx in pybedtools.cbedtools.create_interval_from_list()
    # /tools/anaconda/3-4.4.0/lib/python3.6/site-packages/pybedtools/cbedtools.pyx in pybedtools.cbedtools.isdigit()
    # AttributeError: 'numpy.int64' object has no attribute 'isdigit'
    # SOLUTION: convert everything to strings --> ['chr'+str(x1), str(x2), str(x2)]

    iter_bim = [['chr'+str(x1), str(x2), str(x2)] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)
    
    counter = 1 # just to print status message
    list_df_annot = []
    for name_annotation in sorted(dict_of_beds): # we sort to make output more consistent.
        print("CHR={} | annotation={}, #{}/#{}".format(chromosome, name_annotation, counter, len(dict_of_beds)))
        bed_for_annot = dict_of_beds[name_annotation] # get bed
        # (Pdb)len(bed_for_annot)
        # 66
        #  (Pdb) bed_for_annot.head()
        # chr1    51619935        52185000
        #  chr1   70410488        70871303
        #  chr1   85584164        86243933
        #  chr1   202948059       203355877
        #  chr10  43851792        44270066
        #  chr10  75681524        76110821
        #  chr10  76769912        77191206
        #  chr10  120663598       121138345
        #  chr11  118030300       118469926
        #  chr12  21454715        21871342
        
        annotbed = bimbed.intersect(bed_for_annot, wb=True) # PT NOTE: this finds SNPs in bim file that OVERLAP with the annotation bed (gene)
        # chr22  24008141    24008141    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        # chr22  24008403    24008403    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        # chr22  24008409    24008409    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        # chr22  24008465    24008465    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        # chr22  24008495    24008495    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        # chr22  24008497    24008497    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        # chr22  24008503    24008503    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        # chr22  24008699    24008699    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        # chr22  24008773    24008773    chr22   24008021    24210630    ENSG00000250479 0.03038823367
        
        # annotbed = bed_for_annot.intersect(bimbed, wb=True) # PT NOTE: this finds the positions/intervals in annotation bed (gene) that OVERLAP with the bim file. Only the part of the record intersections occurred
            # *IMPORTANT*: bimbed.intersect(bed_for_annot) and bed_for_annot.intersect(bimbed) DOES NOT return the same positions. However, they do return the same number of 'intersected features'. That is, the returned BedTool object as the same length.
            # bed_for_annot.intersect(bimbed) returns features that span two bp (e.g. start=24008140, end=24008142), whereas bimbed.intersect(bed_for_annot) returns features that span a single bp (start=24008141, end=24008141)
            # use bed_for_annot.intersect(bimbed, wb=True) to understand this behavior better.
        # chr22  24008140    24008142    ENSG00000250479 0.03038823367   chr22   24008141    24008141
        # chr22  24008402    24008404    ENSG00000250479 0.03038823367   chr22   24008403    24008403
        # chr22  24008408    24008410    ENSG00000250479 0.03038823367   chr22   24008409    24008409
        # chr22  24008464    24008466    ENSG00000250479 0.03038823367   chr22   24008465    24008465
        # chr22  24008494    24008496    ENSG00000250479 0.03038823367   chr22   24008495    24008495
        # chr22  24008496    24008498    ENSG00000250479 0.03038823367   chr22   24008497    24008497
        # chr22  24008502    24008504    ENSG00000250479 0.03038823367   chr22   24008503    24008503
        # chr22  24008698    24008700    ENSG00000250479 0.03038823367   chr22   24008699    24008699
        # chr22  24008772    24008774    ENSG00000250479 0.03038823367   chr22   24008773    24008773
        
        ### DOCS .intersect()
        # the intervals reported are NOT the original gene intervals, but rather a refined interval reflecting solely the portion of each original gene interval that overlapped with the SNPs
        # The -wa (write A) and -wb (write B) options allow one to see the original records from the A and B files that overlapped. As such, instead of not only showing you where the intersections occurred, it shows you what intersected.
        # SEE MORE HERE: http://quinlanlab.org/tutorials/bedtools/bedtools.html
        # SEE https://daler.github.io/pybedtools/intersections.html
        # SEE https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.intersect.html#pybedtools.bedtool.BedTool.intersect        
        bp = [x.start for x in annotbed] # PT NOTE: make list of all bp positions for the overlapping SNPs | All features, no matter what the file type, have chrom, start, stop, name, score, and strand attributes.
        annotation_value = [x.fields[7] for x in annotbed] # returns list of strings. Extract the 'score' column. This is column 7 in the 0-based column indexing. *OBS*: x.fields[7] is a string.
        df_annot_overlap_bp = pd.DataFrame({'BP': bp, name_annotation:annotation_value}) # FINUCANE ORIG: df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
        #             BP  blue
        # 0     34605531     1
        # 1     34605604     1
        # 2     34605644     1
        # 3     34605778     1
        # 4     34606634     1
        # 5     34606840     1
        # 6     34607223     1
        df_annot = pd.merge(df_bim, df_annot_overlap_bp, how='left', on='BP') # *IMPORTANT*: how='left' --> resulting data frame will include ALL snps from the bim file.
        # ^ how="left": use only keys from left frame, PRESERVE KEY ORDER
        # (Pdb) df_annot.head()
        #    CHR          SNP        CM       BP  blue
        # 0   21  rs146134162 -0.908263  9412099   NaN
        # 1   21  rs578050168 -0.908090  9412377   NaN
        # 2   21  rs527616997 -0.907297  9413645   NaN
        # 3   21  rs544748596 -0.906578  9414796   NaN
        # 4   21  rs528236937 -0.906500  9414921   NaN
        df_annot = df_annot[[name_annotation]] # get rid of all columns but the name_annotation. Important: return 1 column data frame (and not series, which would loose the column name)
            # df[[name_annotation]] or df.loc[:, [name_annotation]] --> returns dataframe
            # df[name_annotation] or df.loc[:, name_annotation] --> returns series
        df_annot.fillna(0, inplace=True) # SNPs not in df_annot_overlap_bp will have NA values in name_annotation
        # Do data type conversion AFTER .fillna() to avoid problems with NA (float) that cannot be converted to int.
        if args.flag_encode_as_binary_annotation:
            df_annot[name_annotation] = df_annot[name_annotation].astype(int)
            # ALTERNATIVE if you want to handle NANs: df_annot[name_annotation] = pd.to_numeric(df_annot[name_annotation], errors='raise')
        else:
            df_annot[name_annotation] = df_annot[name_annotation].astype(float)
        list_df_annot.append(df_annot)
        if args.flag_annot_file_per_geneset: # write annot file per annotation per chromosome
            file_out_annot = "{}/{}.{}.{}.annot.gz".format(args.out_dir, args.out_prefix, name_annotation, chromosome) # set output filename. ${prefix}.${chr}.annot.gz
            df_annot.to_csv(file_out_annot, sep="\t", index=False, compression="gzip")
        counter += 1
        # if counter == 4: break
    print("CHR={} | Concatenating annotations...".format(chromosome))
    df_annot_combined = pd.concat(list_df_annot, axis='columns') # stack horizontally (there is no joining on indexes, just stacking)
        # *IMPORTANT*: since we did pd.merge(df_bin, df_annot_overlap_bp) with how='left' the know that ALL dfs in list_df_annot have ALL SNPs in df_bim and the order of the SNPs are preserved.
        # ALTERNATIVELY if you don't want 'thin-annot' use this (i.e. adding 'CHR','SNP','CM','BP' columns): df_annot_combined = pd.concat([df_bim]+list_df_annot, axis='columns') # stack horizontally
    
    # print("CHR={} | Calculating standard deviation for annotations...".format(chromosome))
    # df_annot_sd = pd.DataFrame(df_annot_combined.drop(columns=["CHR", "SNP", "CM", "BP"]).std(), columns=["sd"])
    # df_annot_sd.index.name = "annotation"
    # df_annot_sd["n"] = df_annot.shape[1] # number of SNPs in the data frame. This makes it easier to calculate the combined standard deviation across chromosomes later.
    # file_out_annot_combined_sd = "{}/{}.{}.{}.annot_sd".format(args.out_dir, args.out_prefix, "COMBINED_ANNOT", chromosome) 
    # df_annot_sd.to_csv(file_out_annot_combined_sd, sep="\t", index=True)
    ### Output file
    ### annotation      sd      n
    ### antiquewhite3   0.16847050545855485     5
    ### blue1   0.1197907423131066      5
    ### chocolate       0.0     5

    print("CHR={} | Writing annotations...".format(chromosome))
    file_out_annot_combined = get_annot_file_path(chromosome, args) 
    df_annot_combined.to_csv(file_out_annot_combined, sep="\t", index=False, compression="gzip")
    
    # return (annotbed, df_annot_combined)
    return None



###################################### MAIN ######################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_multi_gene_set', type=str, help='a file of gene names, one line per gene.')
    parser.add_argument('--file_gene_coord', type=str, help='a file with columns GENE, CHR, START, and END, where START and END are base pair coordinates of TSS and TES. This file can contain more genes than are in the gene set. We provide ENSG_coord.txt as a default.')
    parser.add_argument('--windowsize', type=int, help='how many base pairs to add around the transcribed region to make the annotation? Finucane uses 100 kb.')
    parser.add_argument('--bimfile_basename', type=str, help='plink bim BASENAME for the dataset you will use to compute LD scores. If argument is "1000G.EUR.QC", then the files "1000G.EUR.QC.1.bim", "1000G.EUR.QC.2.bim", ..., "1000G.EUR.QC.22.bim" will be loaded')
    parser.add_argument('--out_dir', type=str, help='output directory to write annot files. Relative or absolute. Dir be created if it does not exist. ')
    parser.add_argument('--out_prefix', type=str, help='Prefix for output files. Outputfiles will be <out_dir>/<out_prefix>.<name_annotation>.<chromosome>.annot.gz')
    parser.add_argument('--flag_annot_file_per_geneset', action='store_true', help='set flag to write one annot file per gene set per chromosome. Default is to write a combined annot file per chromosome containing all gene sets. NB: the annotation files are always split per chromosome.')
    parser.add_argument('--flag_encode_as_binary_annotation', action='store_true', help='set flag if LDSC annotations should be encoded as binary annotations {0,1}. The default is to use the continuous annotations, which require an appropriate field in the file_multi_gene_set')
    parser.add_argument('--flag_wgcna', action='store_true', help='set flag if file_multi_gene_set input is from WGCNA pipeline')
    parser.add_argument('--flag_mouse', action='store_true', help='set flag if ile_multi_gene_set input has mouse genes (instead of human)')
    parser.add_argument('--n_parallel_jobs', type=int, default=22, help='Number of processes. Default 22 (the number of chromosomes)')

    #parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')


    args = parser.parse_args()

    ### TODO
    # [ ] check that ALL bim files exist. (so pool of workers does not fail)

    # out_prefix = args.out_prefix
    # out_dir = args.out_dir


    ### Make out_dir
    if not os.path.exists(args.out_dir):
        print("Making output dir {}".format(args.out_dir))
        os.makedirs(args.out_dir)


    list_chromosomes = range(1,23) # 1...22
    ### Check for existing annot files
    pool = multiprocessing.Pool(processes=len(list_chromosomes))
    list_chromosomes_to_run = pool.map(functools.partial(check_annot_file, args=args), list_chromosomes) # check_annot_file returns None if annot file exists and is ok
    list_chromosomes_to_run = [x for x in list_chromosomes_to_run if x is not None] # filter away None

    print("========================== CHROMOSOMES TO RUN ====================")
    print("N chromosomes = {}".format(len(list_chromosomes_to_run)))
    print("Chromosomes: {}".format(",".join(map(str, list_chromosomes_to_run))))
    print("=========================================================================")
    if len(list_chromosomes_to_run) == 0:
        print("Quitting...")
        sys.exit()

    ### read coord file
    df_gene_coord = pd.read_csv(args.file_gene_coord, delim_whitespace = True)
    ### read gene list file (a list of gene list files)
    df_multi_gene_set = read_multi_gene_set_file(args)
    ### make beds
    dict_of_beds = multi_gene_sets_to_dict_of_beds(df_multi_gene_set, df_gene_coord, args.windowsize)

    ### estimate number of jobs we can run
    memomory_free_gb = psutil.virtual_memory().available / 1024.0 ** 3 # 1024^3 = Byte to Gigabyte
    n_annotations = df_multi_gene_set["annotation"].nunique()
    ### Memory usage benchmarked Dec 27, 2018
    ### Here we only estimate the memory usage of list_df_annot, which holds n_snps*n_annotations elements
    ### 6000 annotations = 34.8 gb for chr1 (~780k SNPs, float encoded)
    ### 3000 annotations = 17.5 gb for chr1 (~780k SNPs, float encoded)
    ### 1500 annotations = 8.7 gb for chr1 (~780k SNPs, float encoded)
    ### 750 annotations = 4.3 gb for chr1 (~780k SNPs, float encoded)
    ### ---> 5.8 MB per annotation
    memory_needed_per_chromosome_gb = (5.8/1024)*n_annotations # 5.8 MB * n_annotations
    n_processes_with_enough_memory = memomory_free_gb/memory_needed_per_chromosome_gb # this is the amount of processes we have enough free memory for
    n_processes_with_enough_memory_incl_buffer = n_processes_with_enough_memory/2.0 # we reserve each process x2 the amount of memory to account for the doubling of memory usage during pandas/numpy concatenation. READ here about concatenation https://github.com/pandas-dev/pandas/issues/14282#issue-178674951
    print("RESOURCES: memomory_free_gb={} | n_processes_with_enough_memory={} | n_processes_with_enough_memory_incl_buffer={}".format(memomory_free_gb, n_processes_with_enough_memory, n_processes_with_enough_memory_incl_buffer))
    n_parallel_jobs_auto_configured = max(1, min(22 , n_processes_with_enough_memory_incl_buffer))
    n_parallel_jobs_auto_configured

    print("Starting pool with processes={}".format(n_parallel_jobs_auto_configured))
    pool = multiprocessing.Pool(processes=n_parallel_jobs_auto_configured)
    pool.map(functools.partial(make_annot_file_per_chromosome, dict_of_beds=dict_of_beds, args=args), list_chromosomes_to_run)
    # ^ this works for Python 2.7+. (Only Python 3.3+ includes pool.starmap() which makes it easier)
    # ^ Partial creates a new simplified version of a function with part of the arguments fixed to specific values.
    # REF "Python multiprocessing pool.map for multiple arguments": https://stackoverflow.com/a/5443941/6639640
    # pool.map(make_annot_file_per_chromosome_star, itertools.izip(list_chromosomes_to_run, itertools.repeat(dict_of_beds), itertools.repeat(args))) # ALTERNATIVE APPROACH, but requires making a 'make_annot_file_per_chromosome_star()' function


    print("Script is done!")


###################################### Memory usage ######################################

### Code for estimating memory usage of list_df_annot
# n_snps = 780000 # average number of SNPs in .bim files is 450k (150k-780k)
# n_annotations = 750
# x_array = np.arange(n_snps).astype(float)
# df = pd.DataFrame(np.array([x_array]*n_annotations).T)
# mem_usage_gb = df.memory_usage(index=True, deep=True).sum()/(1024.0**3) # deep=If True, introspect the data deeply by interrogating object dtypes for system-level memory consumption, and include it in the returned values.
# mem_usage_gb

