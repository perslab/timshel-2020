

# ======================================================================= #
# =========== PLOT: single gene, facet_wrap over ES metrics ============= #
# ======================================================================= #
### WORKS

# ### 
# df.plot <- df.es
# # df.plot <- df.plot %>% filter(es_metric=="tstat")
# # df.plot <- df.plot %>% filter(es_metric == "mean nES")
# df.plot <- df.plot %>% filter(es_weight > 0)
# p <- ggplot(df.plot, aes(x=es_weight)) + 
#   geom_density() + 
#   scale_x_log10() +
#   facet_wrap(~es_metric, scales = "free")
# p

# ======================================================================= #
# ============================ COLOR density ========================= #
# ======================================================================= #


# =======================  ======================= #
### REF: https://stackoverflow.com/questions/3494593/shading-a-kernel-density-plot-between-two-points
set.seed(1)
draws <- rnorm(100)^2
q75 <- quantile(draws, .75)
q95 <- quantile(draws, .95)
dens <- density(draws)
dd <- with(dens,data.frame(x,y)) # 'dens' contains elements 'x' and 'y'
qplot(x,y,data=dd,geom="line")+
  geom_ribbon(data=subset(dd,x>q75 & x<q95),aes(ymax=y),ymin=0,
              fill="red",colour=NA,alpha=0.5)