

### Load h2 (observed scaled) estimates
file.h2 <- here("results/h2_trait.multi_gwas.csv")
df.h2 <- read_csv(file.h2)

# ======================================================================= #
# ================================ Effect of liability scale h2 barplot =========================== #
# ======================================================================= #
filter.gwas <- c("SCZ_Pardinas2018_liability_scale",
                 "INSOMNIA_Jansen2018_liability_scale",
                 "MS_Patsopoulos2011_liability_scale",
                 "RA_Okada2014_liability_scale",
                 "SCZ_Pardinas2018",
                 "INSOMNIA_Jansen2018",
                 "MS_Patsopoulos2011",
                 "RA_Okada2014")

### Create 'plotting' data frame
df.plot.h2 <- df.h2 %>% filter(gwas %in% filter.gwas)

### Plot
p.h2 <- ggplot(df.plot.h2, aes(x=gwas, y=h2)) + 
  geom_col(fill="gray") +
  labs(y=expression(h[SLDSC]^{2})) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
p.h2