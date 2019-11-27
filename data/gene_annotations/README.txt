
###################################### Uniqueness of ortholog mapping ######################################
timshel@yggdrasil:/nfsdata/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations> zcat gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz | wc -l
16229

ensembl_gene_id=1
timshel@yggdrasil:/nfsdata/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations> zcat gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz | cut -f 1 | sort -u | wc -l
16229

mmusculus_homolog_ensembl_gene=5
timshel@yggdrasil:/nfsdata/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations> zcat gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz | cut -f 5 | sort -u | wc -l
16229