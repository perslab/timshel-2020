
library(tidyverse)

# ======================================================================= #
# =================  TMP cmp z_score vs z_score_two_step  =============== #
# ======================================================================= #

df.z1 <- read_csv("/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain.celltype_expr.z_score.hsapiens_orthologs.csv.gz") %>% column_to_rownames("gene") %>% as.data.frame()
df.z2 <- read_csv("/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain.celltype_expr.z_score_two_step.hsapiens_orthologs.csv.gz") %>% column_to_rownames("gene") %>% as.data.frame()
i <- 5
qplot(df.z1[,i], df.z2[,i]) # cell-type
qplot(as.numeric(df.z1[i,]), as.numeric(df.z2[i,])) # gene


