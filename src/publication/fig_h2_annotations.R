############### SYNOPSIS ###################
# Heritability explained figure


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/publication"))

# ======================================================================= #
# ============================ CREATE DF =============================== #
# ======================================================================= #

df <- read_tsv(
  "idx_sort	name	annotation	gwas	Prop._SNPs	Prop._h2	Prop._h2_std_error	Enrichment
  1	Limb Muscle - mesenchymal stem cell (Height, Loh 2018)	Limb_Muscle.mesenchymal_stem_cell	HEIGHT_UKBB_Loh2018	0.097142921	0.186825909	0.009437455	1.923206618
  2	Liver - hepatocyte (LDL, Teslovich 2010)	Liver.hepatocyte	LIPIDS_LDL_Teslovich2010	0.066339952	0.206528415	0.028341378	3.113183083
  3	Spleen - T cell (Rheumatoid arthritis, Okada 2014)	Spleen.T_cell	RA_Okada2014	0.048971844	0.198452016	0.023819019	4.052369674
  4	Pancreas - Beta cell (T2D, Loh 2018)	Pancreas.type_B_pancreatic_cell	T2D_UKBB_Loh2018	0.15510798	0.222278292	0.014808791	1.433055162
  5	TEGLU23 (BMI, Loh 2018)	TEGLU23	BMI_UKBB_Loh2018	0.182497539	0.241670025	0.006674005	1.324237173
  6	DEINH3 (BMI, Loh 2018)	DEINH3	BMI_UKBB_Loh2018	0.134163963	0.18542722	0.005990666	1.382094086
  7	MEGLU1 (BMI, Loh 2018)	MEGLU1	BMI_UKBB_Loh2018	0.134102478	0.182333335	0.005604643	1.359656713
  8	MEINH2 (BMI, Loh 2018)	MEINH2	BMI_UKBB_Loh2018	0.139807026	0.195322196	0.007040324	1.397084268
  9	DEGLU5 (BMI, Loh 2018)	DEGLU5	BMI_UKBB_Loh2018	0.117483223	0.168099447	0.007425705	1.43083789
  10	MEGLU10 (BMI, Loh 2018)	MEGLU10	BMI_UKBB_Loh2018	0.119019706	0.159239148	0.005135935	1.337922546
  11	TEGLU17 (BMI, Loh 2018)	TEGLU17	BMI_UKBB_Loh2018	0.19148104	0.248760279	0.007483996	1.299137915
  12	MEGLU11 (BMI, Loh 2018)	MEGLU11	BMI_UKBB_Loh2018	0.124795507	0.169261686	0.006306219	1.356312344
  13	TEGLU4 (BMI, Loh 2018)	TEGLU4	BMI_UKBB_Loh2018	0.192775762	0.25788889	0.008282548	1.337766153
  14	DEGLU4 (BMI, Loh 2018)	DEGLU4	BMI_UKBB_Loh2018	0.098469716	0.132322351	0.005001685	1.343787283
  15	TEINH12 (BMI, Loh 2018)	TEINH12	BMI_UKBB_Loh2018	0.13346096	0.179357864	0.006380445	1.3438976"
)

### order
df <- df %>% mutate(name=factor(name, levels=name[order(-idx_sort)]))
df$name %>% levels()

### add percentage
df <- df %>% mutate(Pct._h2=Prop._h2*100)

# ======================================================================= #
# ============================ PLOT =============================== #
# ======================================================================= #


### SIMPLE - with h2% label
ggplot(df, aes(x=name, y=Pct._h2, label=paste0(round(Pct._h2, 1), "%"))) + 
  geom_segment(aes(x=name, xend=name, y=0, yend=Pct._h2), color="grey") +
  geom_point(aes(color=gwas), size=10) + 
  geom_text(size=3) +
  coord_flip() + 
  labs(x="", y=expression("%"~h[2]~explained)) +
  guides(color=F)
ggsave("out.plot.fig_h2.simple_with_labels.pdf", width=10, height=5)


### FULL COMPLEXITY (enrichment + Prop._SNPs)
ggplot(df, aes(x=name, y=Pct._h2)) + 
  geom_point(aes(size=Enrichment, color=gwas)) + 
  geom_segment(aes(x=name, xend=name, y=0, yend=Pct._h2), color="grey") +
  geom_segment(aes(x=name, xend=name, y=0, yend=Prop._SNPs), color="black") +
  coord_flip() + 
  labs(size=expression(h[2]~enrichment), x="", y=expression("%"~h[2]~explained)) +
  guides(color=F)
ggsave("out.plot.fig_h2.complex.pdf", width=10, height=5)



# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


