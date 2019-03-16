### Copy of specificity.index() function from CRAN pSI package.
### Modifications: 
# 1: return both pSI and SI values as a list [no SI = FALSE argument]

### TODO
### 1. Parallelize.
### 2. Make print statements
### 3. Fix issue with if ((!is.na(avg_ipO.i)) & (avg_ipO.i < ctoff)) { # PT NOTE: may cause error 'missing value where TRUE/FALSE needed'
### ---> don't know how to fix it.
# x <- NA
# if ((!is.na(x)) & (x < 5)) {
# print("BLA")
# }



specificity.index.timshel <- function (pSI.in, pSI.in.filter, bts = 50, p_max = 0.1, e_min = 0.3, hist = FALSE) {
  dat_clean <- pSI.in
  dat_clean[pSI.in < e_min] <- NA
  if (!missing(pSI.in.filter)) {
    dat_clean[is.na(pSI.in.filter)] <- NA
    rm(pSI.in.filter)
  }
  lng <- length(pSI.in[1, ])
  datComb <- array(NA, c(length(pSI.in[, 1]), 2 * lng))
  colnames(datComb) <- c(rep(colnames(pSI.in), each = 2))
  psi_columns <- c(seq(from = 2, to = ncol(datComb), by = 2))
  si_columns <- c(seq(from = 1, to = ncol(datComb), by = 2))
  rownames(datComb) <- rownames(pSI.in) # PT note: rownames are used
  for (j in 1:lng) { # j: index cell-types (columns of input)
    notme = c(1:(j - 1), (j + 1):lng)
    if (j == 1) {
      notme = c(2:lng)
    }
    if (j == lng) {
      notme = c(1:(lng - 1))
    }
    TmpMat <- pSI.in[!is.na(dat_clean[, j]), ]
    smps <- length(TmpMat[, 1])
    avg_ip_pre <- array(NA, c(smps, bts))
    for (k in 1:bts) {
      TmpMatSmp <- array(NA, c(smps, lng))
      for (i in 1:lng) {
        TmpMatSmp[, i] <- sample(TmpMat[, i], smps, 
                                 replace = TRUE)
      }
      z <- log2(TmpMatSmp[, j]/TmpMatSmp[, notme])
      rankz <- array(0, dim(z))
      z <- z * -1
      for (i in 1:(lng - 1)) {
        rankz[, i] = rank(z[, i], na.last = "keep")
      }
      avg_ip <- c()
      avg_ip <- .rowMeans(rankz, m = smps, n = ncol(rankz))
      avg_ip_pre[, k] <- avg_ip # PT this is the actual 'return value' for each 'background'/permutation. 
      # No other values in the 'for (k in 1:bts)' loop is used
    } # END 'for (k in 1:bts)' loop
    avg_ip_dist <- sort(avg_ip_pre)
    zO <- array(NA, c(length(dat_clean[, 1]), lng - 1))
    zO <- log2(dat_clean[, j]/pSI.in[, notme])
    rankzO <- array(0, dim(zO))
    zO <- zO * -1
    for (i in 1:(lng - 1)) {
      rankzO[, i] = rank(zO[, i], na.last = "keep")
    }
    avg_ipO <- c()
    avg_ipO <- .rowMeans(rankzO, m = nrow(rankzO), n = ncol(rankzO))
    if (hist == TRUE) {
      fnm = paste("hist_rnks", colnames(pSI.in)[j], ".png", 
                  sep = "_")
      png(filename = fnm, width = 960, height = 960)
      par(mfrow = c(2, 2))
      hist(avg_ip_dist, freq = FALSE, xlim = c(0, smps), 
           main = "Sampled distribution")
      hist(avg_ipO, freq = FALSE, xlim = c(0, smps), main = "Actual distribution")
      dev.off()
    }
    pvals <- c()
    lng2 <- length(avg_ipO)
    ctoff <- avg_ip_dist[(smps * bts) * p_max]
    avg_ip_dist_filt <- avg_ip_dist[which(avg_ip_dist < 
                                            ctoff)]
    for (i in 1:lng2) {
      avg_ipO.i <- avg_ipO[i]
      if ((!is.na(avg_ipO.i)) & (avg_ipO.i < ctoff)) { # PT NOTE: may cause error 'missing value where TRUE/FALSE needed' when p_max=1 and e_min=0
        avg_ip_dist_filt[1] <- avg_ipO.i
        pvals[i] <- ((rank(avg_ip_dist_filt, na.last = "keep")[1])/(smps * 
                                                                      bts))
      }
      else {
        pvals[i] <- NA
      }
    }
    datComb[, ((j * 2) - 1)] <- avg_ipO # SI?
    datComb[, (j * 2)] <- pvals # pSI?
  } # END LOOP 'for (j in 1:lng)'. Index cell-types (columns of input)
  # if (SI == TRUE) {
  #   datComb <- datComb[, si_columns]
  # }
  # else {
  #   datComb <- datComb[, psi_columns]
  # }
  df.SI <- data.frame(datComb[, si_columns])
  df.pSI <- data.frame(datComb[, psi_columns])
  list.res <- list(SI=df.SI, pSI=df.pSI) # named list
  
  return(list.res)
}


 
# library(doParallel)
# library(foreach) # lib is also loaded as part of doParallel, but we just make sure to load it. Needed for "%dopar%"
# registerDoParallel(30)
# time.start <- proc.time()
# iter_control <- df.cem.top %>% select(-cell_type) %>% colnames()
# # iter_control <- iter_control[1:10]
# list.res <- foreach (idx.foreach=1:length(iter_control), # columns
#                      .inorder=TRUE, # .inorder=FALSE can increase performance, but we want things to be in order so we can set the names
#                      .packages=c("tidyverse")) %dopar% {
#                        #list.par_analysis <- foreach (i=foreach_min:foreach_max, .packages=c("plyr")) %dopar% {
#                        time.loop.start <- proc.time()
#                        name_element.foreach <- iter_control[idx.foreach] # an element in iter_control, e.g. "slateblue4"
#                        print(paste("processing #", idx.foreach, "/#", length(iter_control), sep="" ))
#                        print(paste("name_element.foreach is:", name_element.foreach))
#                        df.res.t <- data.frame(cell_type=unique(df.cem.top$cell_type), t_test=NA, pval=NA, df=NA)
#                        for (cell_type_loop in unique(df.cem.top$cell_type)) {
#                          print(cell_type_loop)
#                          x <- df.cem.top %>% filter(cell_type==cell_type_loop) %>% select(UQ(rlang::sym(name_element.foreach)))
#                          y <- df.cem.top %>% filter(cell_type!=cell_type_loop) %>% select(UQ(rlang::sym(name_element.foreach)))
#                          t <- t.test(x,y, alternative="greater",var.equal=T)
#                          df.res.t[df.res.t$cell_type==cell_type_loop, "t_test"] <- t$statistic
#                          df.res.t[df.res.t$cell_type==cell_type_loop, "pval"] <- t$p.value
#                          df.res.t[df.res.t$cell_type==cell_type_loop, "df"] <- t$parameter
#                        }
#                        df.res.t
#                        
#                        # df.x <- df.cem.top %>% 
#                        #   group_by(cell_type) %>% 
#                        #   select(UQ(rlang::sym(name_element.foreach))) %>% # select module
#                        #   summarise(var_c=var(UQ(rlang::sym(name_element.foreach))),
#                        #             mean_c=mean(UQ(rlang::sym(name_element.foreach)))
#                        #             )
#                        # ### "return value"
#                        
#                      }
# ### Setting names of result list
# # foreach(..., .inorder=TRUE, ...): logical flag indicating whether the .combine function requires the task results to be combined in the same order that they were submitted. If the order is not important, then it setting .inorder to FALSE can give improved performance. The default value is TRUE.
# # ---> since the ".inorder" is true per default, the list is combined in the same order that the elements were submitted.
# names(list.res) <- iter_control
# ### Binding
# df.res.foreach <- bind_rows(list.res, .id="module_id")
# df.res.foreach