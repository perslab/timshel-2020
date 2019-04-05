

import sys
import subprocess
import os
import glob

import argparse
import pandas as pd

### Usage
# python run_magma_geneset.py \
# --file_gwas_raw /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw \
# --file_wgcna_module_genelist /projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv \
# --outprefix magma_geneset_BMI_Yengo2018



### SNIPPET WGCNA genelist
# cell_cluster,module,ensembl,hgcn,pkIM
# AP,darkkhaki,ENSMUSG00000041607,Mbp,0.139703362227531
# AP,darkkhaki,ENSMUSG00000031425,Plp1,0.137452357600197
# AP,darkkhaki,ENSMUSG00000036634,Mag,0.137127194286011
# AP,darkkhaki,ENSMUSG00000032517,Mobp,0.133225515039086
# AP,darkkhaki,ENSMUSG00000073680,Tmem88b,0.132201991724902


###################################### FUNCTIONS ######################################
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
    
    
    file_mapping = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
    ### SNIPPET
    # ensembl_gene_id chromosome_name start_position  end_position    mmusculus_homolog_ensembl_gene  mmusculus_homolog_orthology_confidence
    # ENSG00000138593 15      49280673        49338760        ENSMUSG00000035093      1
    # ENSG00000166351 21      14982498        15013906        ENSMUSG00000095294      0
    # ENSG00000168675 18      13217497        13652754        ENSMUSG00000024544      1
    
    df_mapping = pd.read_csv(file_mapping, delim_whitespace=True, usecols=["ensembl_gene_id", "mmusculus_homolog_ensembl_gene"], index_col=1)
    
    ### map ortholog genes
    df["gene_human"] = df["gene"].map(df_mapping["ensembl_gene_id"]) # df["gene"] is mouse genes. df_mapping["ensembl_gene_id"] is human.
    # ^ .map() returns NaN values for genes not mapped. 

    ### make summary of mapping
    file_out_mapping_stats = "{}.ortholog_mapping_stats.txt".format(outprefix)
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


def map_hs_ensembl_to_entrez(df):
    # Entrezgene to ensemb_gene_id mapping file path for MAGMA
    file_mapping_hs_ensembl2entrez = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz"
    df_mapping_ensembl2entrez = pd.read_csv(file_mapping_hs_ensembl2entrez, sep="\t") # shape=(64719, 3)
    df_mapping_ensembl2entrez.drop_duplicates('ensembl_gene_id', keep='first', inplace=True) # shape=(63967, 3)
    df_mapping_ensembl2entrez.drop_duplicates('entrezgene', keep='first', inplace=True) # shape=(25233, 3)
    df_mapping_ensembl2entrez.set_index("ensembl_gene_id", inplace=True) # set index so '.map()' will work
    
    ### map
    df["gene_entrez"] = df["gene"].map(df_mapping_ensembl2entrez["entrezgene"])
    # ^ .map() returns NaN values for genes not mapped. 
    
    ### summarise mapping stats (optional)
    file_out_mapping_stats = "{}.entrez2ensembl_mapping_stats.txt".format(outprefix)
    df_summary = df.groupby("annotation")["gene_entrez"].agg({'n_genes_input': lambda x: len(x),
                                                         'n_genes_output': lambda x: len(x)-sum(pd.isnull(x)), 
                                                        'n_genes_not_mapped' : lambda x: sum(pd.isnull(x)),
                                                        'pct_genes_not_mapped': lambda x: "{:.2f}".format(sum(pd.isnull(x))/float(len(x))*100)})
    # ^ we use pd.isnull() instead of np.isnan() because of the issue described here: https://stackoverflow.com/a/36001191/6639640
    df_summary.sort_values(by=['n_genes_output'], inplace=True)
    # print(df_summary)
    df_summary.sort_values(by=['n_genes_output'], inplace=True) 
    df_summary.to_csv(file_out_mapping_stats, sep="\t")

    ### final processing
    df.dropna(axis=0, inplace=True) # remove non-mapped genes
    df['gene_entrez'] = df['gene_entrez'].astype(int) # convert to int (otherwise it will be float). But conversion cannot be done before NaNs have been removed
    # df['gene'] = df['gene_entrez'] # replace mouse genes with human genes
    # df.drop(['gene_entrez'], axis=1, inplace=True) # remove human gene column
    
    return(df)

###################################### SCRIPT ######################################

parser = argparse.ArgumentParser()
parser.add_argument('--file_gwas_raw', type=str, help='Filepath to MAGMA GWAS raw. E.g. /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw')
parser.add_argument('--file_wgcna_module_genelist', type=str, help='Filepath to WGCNA Pipeline genelist ("*_cluster_module_genes.csv file)". E.g. /projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv')
parser.add_argument('--outprefix', type=str, help='Output prefix. Absolute or relative. E.g. /scratch/test_prefix or just myprefix')
args = parser.parse_args()


## args
file_gwas_raw = args.file_gwas_raw
file_wgcna_module_genelist = args.file_wgcna_module_genelist
outprefix = args.outprefix

### args test
# file_gwas_raw = "/scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw"
# file_wgcna_module_genelist = "/projects/timshel/sc-genetics/sc-genetics/out/out.wgcna/nn_lira_sema/tables/nn_lira_sema_per_brain_area_run1_cell_cluster_module_genes.csv"
# outprefix = "magma_test_run_bmi"



df_multi_gene_set = pd.read_csv(file_wgcna_module_genelist, usecols=["cell_cluster", "module", "ensembl"]) 
# ^ 'module' = first column; will later be renamed to 'annotation'
# ^ 'ensembl' = second column; will later be renamed to 'gene'
df_multi_gene_set["cell_cluster"] = df_multi_gene_set["cell_cluster"].str.replace(" ", ":WS:") # converting whitespace to a temporary whitespace placeholder because MAGMA uses whitespace delimted files. Relevant for MACA cell-type names. 'Aorta_endothelial cell-antiquewhite3' to Aorta_endothelial:WS:cell-antiquewhite3
df_multi_gene_set["annotation"] = df_multi_gene_set["cell_cluster"] + "-" + df_multi_gene_set["module"] # create annotation column. E.g "CD47+ B cell - darkkhaki". Note that 'module color name' will NEVER contain hyphens
df_multi_gene_set.rename(columns={'ensembl': 'gene'}, inplace=True) # df.rename(columns={'oldName1': 'newName1', 'oldName2': 'newName2'}, inplace=True)
### map to human
df_multi_gene_set_hs = map_ensembl_genes_mouse_to_human(df_multi_gene_set)
df_multi_gene_set_hs.head()
### map ensembl to entrez
df_multi_gene_set_hs_entrez = map_hs_ensembl_to_entrez(df_multi_gene_set)
df_multi_gene_set_hs_entrez.head()
### write output file
fileout_magma_geneset = "{}.set_annot.txt".format(outprefix)
df_multi_gene_set_hs_entrez[["annotation", "gene_entrez"]].to_csv(fileout_magma_geneset, sep="\t", header=False, index=False) # no header, no index.

### run MAGMA
# magma --gene-results <GWAS>.genes.raw --set-annot <INPUT> --out <GWAS.MODULE>
cmd = "/tools/magma/1.06/magma --gene-results {file_gwas_genes_raw} --set-annot {file_geneset} set-col=1 gene-col=2 --out {fileout_prefix}".format(file_gwas_genes_raw=file_gwas_raw, file_geneset=fileout_magma_geneset, fileout_prefix=outprefix)
print(cmd)

subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)


print("SCRIPT DONE")





