####################################################################################################

# The "official" tutorial is here:
# https://github.com/bulik/ldsc/wiki

# plink_path <- "/Volumes/Seagate_Backup_Plus_Drive/bioinformatics/plink_mac/"
# ldsc_path <- "/Volumes/Seagate_Backup_Plus_Drive/bioinformatics/ldsc/"
# workdir <- "~/Dropbox/projects/new_projects/ldsc/IGSS/tutorial"

# plink_path <- "/usr/local/bin/"
# ldsc_path <- "~/ldsc/"
# workdir <- "~/tutorial"

# ======================= SETUP ======================= #


ldsc_path <- "/projects/timshel/sc-genetics/ldsc/ldsc/"
workdir <- "/projects/timshel/sc-genetics/ldsc/src/tutorial-james/"
python_executive <- "/tools/anaconda/2-4.4.0/bin/python" # OBS: hack: trailing white space.

gwas_data_path <- "/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats_raw/JAMESLEE_traits"
ldsc_data_path <- "/projects/timshel/sc-genetics/ldsc/data"

setwd(workdir)

# ======================= CREATING LD SCORES ======================= #

### PRELIMINARY ###

system(paste(python_executive, " ", ldsc_path, "ldsc.py -h", sep=""))
system(paste(python_executive, " ", ldsc_path, "munge_sumstats.py -h", sep=""))

### CREATING LD SCORES ###

system(sprintf("zcat %s | head | column -t", paste(ldsc_data_path, "/eur_ref_ld_chr/22.l2.ldscore.gz", sep="")))
# The column CM refers to centimorgan, which is the distance between CHR positions for which the
# expected number of intervening chromosomal crossovers in a single generation is 0.01. In
# humans this is about 1 Mb on average.

# The precomputed LD Scores (based on 1000G) available from
# https://data.broadinstitute.org/alkesgroup/LDSCORE/
# use a window of 1 CM around each HM3 SNP and exclude singletons (MAF \approx 0.0013). These
# should be adequate for most uses. Incidentally, this script assumes that all of the files
# accessible at the URL above have been placed in their subdirectories in "files," itself
# a subdirectory of the LDSC path.

# There might be situations where you want to compute your own LD Scores. In such a case you will
# probably download the 1000G VCF files from
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
# and convert them to PLINK format (plink --vcf). Select individuals from a given population 
# (e.g., East Asians) using a key such as the one you can find by googling
# "20101214_1000genomes_samples".

# WARNING: The PLINK files that you create in this way will not have the CM information! You have
# to get this from somewhere else, such as
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# In this tutorial we will just compute LD Scores of SNPs on chromosome 22, using only HM3 SNPs.
# Usually you will want to use 1000G, but just HM3 SNPs can be useful for determining the
# regression weights.


command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --bfile ", ldsc_data_path, "/1kg_eur/22",
                 " --l2",
                 " --ld-wind-cm 1",
                 " --out 22",
                 sep="")
# If you do not have a true CM column, you can use --ld-wind-kb or --ld-wind-snp instead.                 
system(command)                 
# There should be a correlation of ~0.3 btw MAF and LDSC.

# A rough sanity check is to calculate the correlations between your new LD Scores and old ones
# that you already have. Even if they are calculated differently, there should be some correlation,
# especially if the population is the same.

new_ldsc <- read.table("22.l2.ldscore.gz", header=TRUE, colClasses=c("integer", "character",
                       "integer", "numeric"))
head(new_ldsc)

old_ldsc <- read.table(paste(ldsc_data_path, "/eur_w_ld_chr/22.l2.ldscore.gz", sep=""),
                       header=TRUE, colClasses=c("integer", "character", "integer", "numeric",
                       "numeric", "numeric"))                       
head(old_ldsc)

both_ldsc <- merge(old_ldsc, new_ldsc, by="SNP")
cor.test(both_ldsc$L2.x, both_ldsc$L2.y)  # ---> 0.9998879 

# system("rm -f 22.*")




# ======================= EA2 - .meta ======================= #
# ### ESTIMATING INTERCEPTS AND CORRECTING STANDARD ERRORS ###
# 
# system(sprintf("zcat %s/EA2_excl_23andMe_noGC.meta.gz | head", gwas_data_path))
# 
# command <- paste(python_executive, " ", ldsc_path, "munge_sumstats.py",
#                  " --sumstats ", gwas_data_path, "/EA2_excl_23andMe_noGC.meta.gz",
#                  " --out EA2_excl_23andMe_noGC",
#                  " --merge-alleles ", ldsc_data_path, "/w_hm3.snplist",
#                  " --ignore Direction,MarkerName",
#                  sep="")
# # GWAS summary stats will generally not be in the format required by LDSC regression. The
# # munge_sumstats command reformats the summary stats and performs a number of useful checks.
# # The --merge-alleles option ensures that your SNPs have the alleles listed in HM3.
# 
# system(command)             
# 
# # If munge_sumstats does not work for some reason, you can always manually reformat the GWAS
# # summary stats. Let's take a look at what the LDSC regression format looks like.
# 
# system("gzcat EA2_excl_23andMe_noGC.sumstats.gz | head")
# 
# # Now let's look at the "heritability" and the intercept.
# 
# command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
#                  " --h2 EA2_excl_23andMe_noGC.sumstats.gz",
#                  " --ref-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
#                  " --w-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
#                  " --out EA2_excl_23andMe_noGC",
#                  sep="")
# # The --ref-ld-chr flag specifies which LD Scores to use as the independent variable in the
# # regression. The w-ld-chr flag specifies which LD Scores to use in the weighting of each data
# # point. The LDSC papers and the official tutorial do not give internally consistent 
# # recommendations regarding the LD Scores for weighting. Here I just use the same LD Scores as
# # independent variables and as weights because the choice of weights does not strongly affect
# # the results; this is because LD Scores specified in different ways will be strongly
# # correlated.
# 
# system(command)                             
# 
# # The intercept is 1.11, which might suggest a fair degree of confounding. But this
# # meta-analysis includes a sample of >40K contributed by deCODE, where many individuals are
# # closely related. Thus the effective sample size is not as large as it seems.
# 
# # We can use the intercept to correct for both confounding and the non-independence of the
# # related individuals and thereby achieve a correct Type 1 error rate among null SNPs. 
# 
# system("gzcat EA2_excl_23andMe_noGC.meta.gz | head -n1")
# EA_meta <- read.delim("EA2_excl_23andMe_noGC.meta.gz", colClasses=c(
# "character", # rsID
# "character", # MarkerName
# "NULL", # Allele1
# "NULL", # Allele2
# "NULL", # Freq1
# "NULL", # FreqSE
# "NULL", # MinFreq
# "NULL", # MaxFreq
# "NULL", # Weight
# "numeric", # Zscore
# "numeric", # P-value
# "NULL", # Direction
# "NULL", # HetISq
# "NULL", # HetChiSq
# "NULL", # HetDf
# "NULL" # HetPVal
# ))
# 
# EA_meta <- EA_meta[1:1000000,]
# 
# ChiSq <- EA_meta$Zscore^2
# ChiSqCorrect <- ChiSq/1.11
# ZscoreCorrect <- sqrt(ChiSqCorrect)*sign(EA_meta$Zscore)
# indices <- sample(1:nrow(EA_meta), 1E5)
# plot(EA_meta[indices,"Zscore"], ZscoreCorrect[indices])
# abline(a=0, b=1, lwd=3, col="red")
# 
# PvalCorrect <- pnorm(abs(ZscoreCorrect), lower.tail=FALSE)*2
# nrow(subset(EA_meta, P.value<5E-8))
# length(PvalCorrect[PvalCorrect<5E-8])

# ======================= EA2 ======================= #

command <- paste(python_executive, " ", ldsc_path, "munge_sumstats.py",
                 " --sumstats ", gwas_data_path, "/EduYears_Main.txt.gz",
                 " --out EA2_excl_23andMe_noGC",
                 " --merge-alleles ", ldsc_data_path, "/w_hm3.snplist",
                 " --N 217568", # PT calculated
                 sep="")
system(command)

# If munge_sumstats does not work for some reason, you can always manually reformat the GWAS summary stats. 
# Let's take a look at what the LDSC regression format looks like.
# SNP     A1      A2      N       Z
# rs28576697      T       C       236480.000      -1.344
# rs1110052       T       G       232068.000      -2.106
# rs7523549
# rs3748592       A       G       245800.000      -1.473

# Now let's look at the "heritability" and the intercept.

command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --h2 EA2_excl_23andMe_noGC.sumstats.gz",
                 " --ref-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --w-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --out EA2_excl_23andMe_noGC",
                 sep="")

# ======================= IQ ======================= #
# Let's look at IQ ("cognitive performance"), infant head circumference, and bipolar disorder.

command <- paste(python_executive, " ", ldsc_path, "munge_sumstats.py",
                 " --sumstats ", gwas_data_path, "/CHIC_Summary_Benyamin2014.txt.gz",
                 " --out IQ",
                 " --merge-alleles ", ldsc_data_path, "/w_hm3.snplist",
                 " --N 17989",
                 " --signed-sumstats EFFECT_A1,0",
                 sep="")
# munge_sumstats was not able to determine the column giving the effect of the reference
# allele without the --signed-sumstats flag. Usually, if a column is missing or munge_sumstats
# cannot find it, the error message will tell you what went wrong.                 
                 
system(command)
command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --h2 IQ.sumstats.gz",
                 " --ref-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --w-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --out IQ", sep="")  
system(command)

# ======================= infant head circumference ======================= #
command <- paste(python_executive, " ", ldsc_path, "munge_sumstats.py",
                 " --sumstats ", gwas_data_path, "/EGG_HC_DISCOVERY.txt.gz",
                 " --out HC",
                 " --merge-alleles ", ldsc_data_path, "/w_hm3.snplist",
                 " --N-col TotalSampleSize",
                 sep="")
system(command)
command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --h2 HC.sumstats.gz",
                 " --ref-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --w-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --out HC", sep="")  
system(command)
# ======================= bipolar ======================= #
command <- paste(python_executive, " ", ldsc_path, "munge_sumstats.py",
                 " --sumstats ", gwas_data_path, "/pgc.cross.BIP11.2013-05.txt.gz",
                 " --out BIP",
                 " --merge-alleles ", ldsc_data_path, "/w_hm3.snplist",
                 " --N-cas 6990",
                 " --N-con 4820",
                 sep="")
system(command)
command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --h2 BIP.sumstats.gz",
                 " --ref-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --w-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --out BIP",
                 sep="")  
system(command)

# ======================= GENETIC CORRELATIONS ======================= #
### ESTIMATING GENETIC CORRELATIONS ###

command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --rg EA2_excl_23andMe_noGC.sumstats.gz,IQ.sumstats.gz,",
                 "HC.sumstats.gz,BIP.sumstats.gz",
                 " --ref-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --w-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --out gencorr_EA",
                 sep="")
# This computes the genetic correlation btw the first trait and each of the subsequent traits.                 
                                                                 
system(command)                 

command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --rg IQ.sumstats.gz,HC.sumstats.gz,BIP.sumstats.gz",
                 " --ref-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --w-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --out gencorr_IQ",
                 sep="")
system(command)                 

command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --rg HC.sumstats.gz,BIP.sumstats.gz",
                 " --ref-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --w-ld-chr ", ldsc_data_path, "/eur_ref_ld_chr/",
                 " --out gencorr_HC_BIP",
                 sep="")
system(command) 


# ======================= FUNCTIONAL PARTITION OF HERITABILITY ======================= #


command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --h2 EA2_excl_23andMe_noGC.sumstats.gz",
                 " --ref-ld-chr ", ldsc_data_path, "/1000G_EUR_Phase3_baseline/baseline.",
                 " --w-ld-chr ", ldsc_data_path, "/weights_hm3_no_hla/weights.",
                 " --overlap-annot",
                 " --frqfile-chr ", ldsc_data_path, "/1000G_Phase3_frq/1000G.EUR.QC.",
                 " --print-coefficients",
                 " --out EA2_excl_23andMe_noGC_baseline",
                 sep="")
# The official tutorial recommends using the weights given above. Probably the only important
# consideration here is that the weights represent *total* (non-stratified) LD Scores because
# non-independence and heteroskedasticity of the chi-square statsistics are due to the total
# LD Scores.

# The flag --overlap-annot indicates that annotations are not disjoint. E.g., a SNP can be
# annotated as both coding and evolutionarily conserved.

# The flag --frqfile_chr is used to figure out which SNPs have MAF > 0.05. By default LD Score
# regression uses only these SNPs to extrapolate h^2 because most HM3 SNPs are common,
# whereas there are many more rare than common SNPs. Common and rare SNPs probably do not have
# the same h^2 per SNP, and thus avoiding the extrapolation makes sense. The official tutorial has
# an explanation of why the .frq files must be explicitly flagged when partitioning heritability
# that I find admittedly difficult to parse.

# By default the coefficients of the annotations are not included in the output. The flag
# --print-coefficients overrides this.
                 
system(command)    

baseline <- read.table("EA2_excl_23andMe_noGC_baseline.results", header=TRUE)
head(baseline)

baseline[order(-baseline[,5]),c("Category", "Prop._SNPs", "Enrichment", "Enrichment_p")]
# Different traits tend to yield similar results. For example, just about all GWAS traits
# shows > 10-fold enrichment of evolutionarily conserved regions.

# Different traits tend to yield dramatically different results when we look at tissue-specific
# histone marks.

tissues <- read.table(paste(ldsc_data_path, "/1000G_Phase3_cell_type_groups/names", sep=""),
                    header=TRUE)

for(i in 1:10){
	print(i)
	tissue <- tissues[i,2]
	command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
	                 " --h2 EA2_excl_23andMe_noGC.sumstats.gz",
	                 " --ref-ld-chr ", ldsc_data_path, "/1000G_Phase3_cell_type_groups/",
	                 "cell_type_group.", i, ".,",
	                 ldsc_data_path, "/1000G_EUR_Phase3_baseline/baseline.",
	                 " --w-ld-chr ", ldsc_data_path, "/weights_hm3_no_hla/weights.",
	                 " --overlap-annot",
	                 " --frqfile-chr ", ldsc_data_path, "/1000G_Phase3_frq/1000G.EUR.QC.",
	                 " --print-coefficients",
	                 " --out EA2_excl_23andMe_noGC_", tissue,
	                 sep="")
	system(command) }
# The loop above takes more than 30 min on my desktop. It can be parallelized, but the options
# for doing this vary across operating systems. We'll just look at the precomputed results.

system("head EA2_excl_23andMe_noGC_Adrenal_Pancreas.bed.results")

command <- paste("head -n1 EA2_excl_23andMe_noGC_Adrenal_Pancreas.bed.results >",
                   " EA2_excl_23andMe_noGC_tissue.results")
system(command)
for(i in 1:10){
	tissue <- tissues[i,2]
	command <- paste("head -n2 EA2_excl_23andMe_noGC_", tissue, ".results |",
	                 " tail -1 >> EA2_excl_23andMe_noGC_tissue.results",
	                 sep="")
    system(command) }	
    
tissue_results <- read.table("EA2_excl_23andMe_noGC_tissue.results", header=TRUE)
tissue_results[,1] <- tissues[,2]
barx <- barplot(tissue_results$Coefficient, ylim=c(-5e-08, 5e-08),
        ylab="Stratified LDSC regression tissue coefficients", col=c("red", "orange", "yellow",
        "green", "blue", "purple", "brown", "pink", "gray", "seagreen"))
legend("topright", legend=tissue_results$Category, ncol=2, cex=1, fill=c("red", "orange",
       "yellow", "green", "blue", "purple", "brown", "pink", "gray", "seagreen"))        
abline(h=0)

error_bar <- function(x, y, upper, lower=upper, length=0.1, ...){
	arrows(x, y+upper, x, y-lower, angle=90, code=3, length=length, ...) }
	
error_bar(barx, tissue_results$Coefficient, 2*tissue_results$Coefficient_std_error)

# ======================= MAKING YOUR OWN ANNOTATION ======================= #
# What if I want to test my own annotation? I will demonstrate with DNase I hypersensitivity in
# the fetal brain (DS15453), the top fgwas result when applied to years of education (Okbay et al.,
# 2016). fgwas and LD Score regression are somewhat similar in that both tell us whether an
# SNPs with a particular annotation are more strongly associated with the trait.

# Use the CNS annotation file as a template.

cc <- c("character", "integer", "character", rep("NULL", 451))
cc[c(90,382)] <- "integer"
# fgwas makes use of 451 annotations. It is much quicker to only read in the annotations of
# interest. 

setwd("new_annot")

for(i in 1:22){
	print(i)
	CNS <- read.table(paste(ldsc_data_path, "/1000G_Phase3_cell_type_groups/",
	                  "cell_type_group.3.", i, ".annot.gz", sep=""), header=TRUE,
	                  colClasses=c("integer", "numeric", "character", "numeric", "numeric"))
	fgwas <- read.table(paste("/Volumes/Seagate_Backup_Plus_Drive/bioinformatics/",
	                    "fgwas-0.3.6/annotations/chr", i,
	                    ".annot.wdist.wcoding.gz", sep=""), header=TRUE, colClasses=cc)
    colnames(fgwas)[3] <- "SNP"
    row_order <- 1:nrow(CNS)
    CNS <- data.frame(CNS, row_order)
    
    DNase_fBrain <- merge(CNS[,c(1:4,6)], fgwas[,3:4], by="SNP", all.x=TRUE)
    DNase_fBrain <- DNase_fBrain[ order(DNase_fBrain[,5]), ]
    DNase_fBrain <- DNase_fBrain[,-5]
    colnames(DNase_fBrain)[5] <- "ANNOT"
    DNase_fBrain$ANNOT[is.na(DNase_fBrain$ANNOT)] <- 0
    
    DNase_NPC <- merge(CNS[,c(1:4,6)], fgwas[,c(3,5)], by="SNP", all.x=TRUE)
    DNase_NPC <- DNase_NPC[ order(DNase_NPC[,5]), ]
    DNase_NPC <- DNase_NPC[,-5]
    colnames(DNase_NPC)[5] <- "ANNOT"
    DNase_NPC$ANNOT[is.na(DNase_NPC$ANNOT)] <- 0    

	write.table(DNase_fBrain, paste("DNase_fBrain.", i, ".annot", sep=""), quote=FALSE,
	            sep="\t", row.names=FALSE)
    write.table(DNase_NPC, paste("DNase_NPC.", i, ".annot", sep=""), quote=FALSE,
	            sep="\t", row.names=FALSE) }            
	            
system("gzip DNase_*")
# The code above replaces the final column of the CNS file (the binary indicator of whether
# the SNP bears the annotation) with the desired one while leaving everything else (e.g., the
# order of the SNPs) the same. The fgwas annotation files are very large and not included
# in the tutorial; the code is thus merely illustrative.


# ======================= CALCULATING STRATIFIED LD-scores for NEW annotation ======================= #
for(i in 1:22){
	print(i)
	command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
	                 " --l2",
	                 " --bfile ", ldsc_data_path, "/1000G_EUR_Phase3_plink/1000G.EUR.QC.", i,
	                 " --ld-wind-cm 1",	            
	                 " --annot DNase_fBrain.", i, ".annot.gz",
	                 " --out DNase_fBrain.", i,
	                 " --print-snps ", ldsc_data_path, "/1000G_EUR_Phase3_baseline/print_snps.txt",
	                 sep="")
	system(command)
	
	command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
	                 " --l2",
	                 " --bfile ", ldsc_data_path, "/1000G_EUR_Phase3_plink/1000G.EUR.QC.", i,
	                 " --ld-wind-cm 1",	            
	                 " --annot DNase_NPC.", i, ".annot.gz",
	                 " --out DNase_NPC.", i,
	                 " --print-snps ", ldsc_data_path, "/1000G_EUR_Phase3_baseline/print_snps.txt",
	                 sep="")
	system(command) }	                 
# The code above creates the new stratified LD Scores. The --print-snps flag ensures that only
# the LD Scores of SNPs present in the baseline files are computed; the SNPs in the baseline
# annotation files and the new annotation file must match. This code takes about 2 hours to
# run on my desktop. Let's just work on the precomputed results.


# ======================= Running stratified LD Score regression ======================= #
setwd(workdir)

command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --h2 EA2_excl_23andMe_noGC.sumstats.gz",
                 " --ref-ld-chr new_annot/DNase_fBrain.,",
                 ldsc_data_path, "/1000G_EUR_Phase3_baseline/baseline.",
                 " --w-ld-chr ", ldsc_data_path, "/weights_hm3_no_hla/weights.",
                 " --overlap-annot",
                 " --frqfile-chr ", ldsc_data_path, "/1000G_Phase3_frq/1000G.EUR.QC.",
                 " --print-coefficients",
                 " --out EA2_excl_23andMe_noGC_DNase_fBrain",
                 sep="")
system(command)

command <- paste(python_executive, " ", ldsc_path, "ldsc.py",
                 " --h2 EA2_excl_23andMe_noGC.sumstats.gz",
                 " --ref-ld-chr new_annot/DNase_NPC.,",
                 ldsc_data_path, "/1000G_EUR_Phase3_baseline/baseline.",
                 " --w-ld-chr ", ldsc_data_path, "/weights_hm3_no_hla/weights.",
                 " --overlap-annot",
                 " --frqfile-chr ", ldsc_data_path, "/1000G_Phase3_frq/1000G.EUR.QC.",
                 " --print-coefficients",
                 " --out EA2_excl_23andMe_noGC_DNase_NPC",
                 sep="")
system(command)