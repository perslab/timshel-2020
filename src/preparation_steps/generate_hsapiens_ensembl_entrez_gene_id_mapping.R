############### SYNOPSIS ###################
# Get biomaRt annotation

### DESCRIPTION
# ...

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:
# Ensembl biomaRt guide: http://www.ensembl.org/info/data/biomart/biomart_r_package.html

############################################

### GET MAPPING: ENTREZ --> ENSEMBL
ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
df.mapping <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'entrezgene'), mart = ensembl) # data frame
df.mapping %>% write_delim("gene_id_mapping.hsapiens.ensembl_entrez.txt.gz", delim="\t")


### Our mapping is not unique, but we don't care too much about that at this point.
sum(duplicated(df.mapping$ensembl_gene_id)) # 662 ensembl IDs are duplicated
sum(duplicated(na.omit(df.mapping$entrezgene))) # 3774 entrez IDs are duplicated

df.mapping$ensembl_gene_id[duplicated(df.mapping$ensembl_gene_id)]
df.mapping %>% filter(ensembl_gene_id == "ENSG00000228340")
# Looked up one  gene and it turns out to be some (uncharacterized) microRNA
# > df.mapping %>% filter(ensembl_gene_id == "ENSG00000228340")
# ensembl_gene_id entrezgene
# 1 ENSG00000228340     284757
# 2 ENSG00000228340     729296
