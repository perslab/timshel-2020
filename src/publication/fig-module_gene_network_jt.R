############### SYNOPSIS ###################
### AIM: create a network gene module plot
# JON VERSION

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

library(igraph)
library(WGCNA)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/publication"))

# ======================================================================= #
# ================================ FUNCTION ================================ #
# ======================================================================= #


plotnetwork <- function(module, 
                        moduleColors, 
                        datExpr, 
                        kMEs, 
                        maxsize=40, 
                        path=NULL,
                        returnAdjMat = F) {
  #' @argument module: the gene module to plot, character of length 1
  #' @argument moduleColors: vector of module labels with gene names, character
  #' @argument datExpr: cell * gene matrix, gene colnames must match moduleColors names
  #' @argument kMEs: output of WGCNA::signedKME, numeric data.frame with gene rownames matching moduleColors names and datExpr colnames
  #' @argument maxsize: max number of genes to display, integer, defaults to 40
  #' @argument path: path to which to write plot as pdf. Default is to send plot to default device.
  #' @argument returnAdjMat: return adjacency matrix? Boolean, defaults to F
  #' @value NULL
  #' @usage plotnetwork(module="lavenderblush", moduleColors = moduleColors, datExpr = datExpr, kMEs=kMEs, path=paste0(dirPlots, prefixOut, "_", module, "_networkplot.pdf"), returnAdjMat = F)
  
  # Get the idx and kME values of top module genes ordered by decreasing kME
  gene_idx <- order(kMEs[,module], decreasing = T)[1:maxsize]
  names(gene_idx) <- rownames(kMEs)[gene_idx]
  mod_kMEs = kMEs[gene_idx, module]
  
  # make non-negative correlation matrix -> adjacency matrix
  if (class(all.equal(colnames(datExpr), rownames(kMEs)))=="character") datExpr<-datExpr[,order(rownames(kMEs))]
  adjMat = cor(datExpr[,gene_idx])
  
  # only keep positive correlations
  adjMat[adjMat<0] <- 0
  
  # make graph
  g1 <- graph.adjacency(as.matrix(adjMat),
                        mode="undirected",
                        weighted=T,
                        diag=FALSE)
  edgecolors = numbers2colors(E(g1)$weight, colors = redWhiteGreen(100, gamma=4), signed=T, centered=T, lim=c(-1,1))
  
  #plot
  if (!is.null(path)) pdf(path,useDingbats = F, width=12,height=12)
  plot.igraph(g1, 
              vertex.label = names(gene_idx),
              vertex.label.dist=0, 
              edge.width=0.25,
              vertex.size=mod_kMEs*30, 
              vertex.frame.color="black",
              vertex.label.color="black",
              vertex.color = "white",#modColors[[1]][gene_idx],
              vertex.label.cex=1,
              layout=layout.fruchterman.reingold(g1),
              edge.color="green")
  if(!is.null(path)) invisible(dev.off())
  
  if (returnAdjMat) return(adjMat)
  
}

# ======================================================================= #
# ================================ MAIN ================================ #
# ======================================================================= #


dirFiles = "/projects/jonatan/applied/18-mousebrain_7/tables/networkplot_args/"
dirPlots ="/projects/jonatan/applied/18-mousebrain_7/plots/"
prefixOut ="mb_7.2_190408"

celltype <- "MEINH2" 
#celltype <- "TEGLU23"

module <-  "lavenderblush" 
#module <- "lightpink3"

moduleColors = readRDS(paste0(dirFiles, prefixOut,"_", module, "_", celltype, "_moduleColors.RDS.gz"))
datExpr = read.csv(paste0(dirFiles, prefixOut,"_", module, "_", celltype, "_datExprScaled.csv.gz"), quote="", stringsAsFactors = F, row.names = 1)
kMEs =read.csv(paste0(dirFiles, prefixOut,"_", module, "_", celltype, "_kMEs.csv.gz"), quote="", stringsAsFactors = F, row.names=1)


plotnetwork(module=module, 
            moduleColors = moduleColors, 
            datExpr = datExpr, 
            maxsize = 40,
            kMEs = kMEs,
            path="figs/fig_module_gene_networkplot_jt.pdf",
            returnAdjMat = F)


# TOOD: Figure gene network. ----> Color genes by MAGMA p-val.