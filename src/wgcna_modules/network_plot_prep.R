############ NETWORKPLOTPREP FUNCTION ##############
require(WGCNA)
networkplotPrep <- function(dirWGCNARObjects,
                            prefix_data, 
                            prefix_run,
                            nGenesPlot= 40,
                            module) {
  #' @argument dirWGCNAObjects: path to WGCNA ROBjects subdirectory
  #' @argument prefix_data: data prefix, character
  #' @argument prefix_run: run prefix, character
  #' @argument nGenesPlot: number of genes to plot, defaults to 40, integer
  #' @argument module: the gene module to plot, character
  
  #' @value list with two components:
  #    vec_kMs, vector of module kMs, amed numeric vector
  #    mat_adj, gene-gene adjacency matrix, matrix
  #' @usage list_networkPlotArgs <- networkplotPrep(dirWGCNARObjects = "~/WGCNAanalysis/RObjects/", prefix_data="sc", prefix_run="1", module="yellowgreen")
  
  # Load WGCNA run image
  path_image <- dir(path=dirWGCNARObjects, pattern=paste0(prefix_data, "_", prefix_run, "_final_session_image.RData.gz"), full.names = T)
  newEnv <- new.env()
  load(file=path_image, envir=newEnv)
  
  # get vector of module labels with gene names, character
  df_geneModule <- newEnv$df_meta_module_genes
  colname_kMs <- paste0("p", newEnv$params_run$fuzzyModMembership) #grep("^pk",colnames(df_geneModule),value = T)
  
  # Get ordered, named vec_kMEs
  idx_modGenes <- which(df_geneModule$module==module)
  idx_modGenesOrder <- idx_modGenes[order(df_geneModule[[colname_kMs]][idx_modGenes], decreasing = T)]
  vec_kMs <- df_geneModule[[colname_kMs]][idx_modGenesOrder]
  
  colGenes <- if (newEnv$params_run$map_genes_to_ensembl) "ensembl" else "hgnc"
  names(vec_kMs) <- df_geneModule[[colGenes]][idx_modGenesOrder]
  vec_kMs <- vec_kMs[1:min(length(vec_kMs), nGenesPlot)]
  
  # Get datExpr 
  cellClust <- unique(df_geneModule$cell_cluster[df_geneModule$module==module])
  datExpr <- newEnv$list_datExpr_meta[[cellClust]]
  
  # compute correlation adjacency matrix
  mat_adj = cor(datExpr[,names(vec_kMs)])
  
  # only keep positive correlations
  mat_adj[mat_adj<0] <- 0
  
  return(list("vec_kMs"=vec_kMs, "mat_adj"=mat_adj))
}

############ PLOTGENENETWORK FUNCTION ##############

require(igraph)
require(WGCNA)

plotGeneNetwork <- function(vec_kMs, mat_adj) {
  #' @argument vec_kMs: vector of module kMs, as returned by networkplotPrep, named numeric vector
  #' @argument mat_adj: gene-gene adjacency matrix, as returned by networkplotPrep, matrix
  
  #' @value returns NULL invisibly
  #' @usage 
  #'  pdf(paste0(dirPlots, prefixOut, "_", module, "_networkplot.pdf"),useDingbats = F, width=12,height=12)
  #'  plotGeneNetwork(vec_kMs=list_networkPlotArgs$vec_kMs, mat_adj=list_networkPlotArgs$mat_adj)
  #'  dev.off()
  
  g1 <- graph.adjacency(as.matrix(mat_adj),
                        mode="undirected",
                        weighted=T,
                        diag=FALSE)
  edgecolors = numbers2colors(E(g1)$weight, colors = redWhiteGreen(100, gamma=4), signed=T, centered=T, lim=c(-1,1))
  
  plot.igraph(g1, 
              vertex.label = names(vec_kMs),
              vertex.label.dist=0, 
              edge.width=0.25,
              vertex.size=vec_kMs*30, 
              vertex.frame.color="black",
              vertex.label.color="black",
              vertex.color = "white",#modColors[[1]][gene_idx],
              vertex.label.cex=1,
              layout=layout.fruchterman.reingold(g1),
              edge.color="green")
  #if(!is.null(path)) invisible(dev.off())
}

############# TEST RUN ############
dir_plots=  "/projects/jonatan/applied/18-mousebrain_7/plots/"
module="lavenderblush"
prefix_run="run1"
prefix_data="Neurons_sub_ClusterName_7.2"
dirWGCNARObjects ="/projects/jonatan/applied/18-mousebrain_7/RObjects/"

list_networkPlotArgs <- networkplotPrep(dirWGCNARObjects =dirWGCNARObjects, 
                                        prefix_data=prefix_data, 
                                        prefix_run = prefix_run, 
                                        nGenesPlot = 50,
                                        module=module)

pdf(paste0(dir_plots, prefix_data, "_", prefix_run, "_", module, "_networkplot_test.pdf"),useDingbats = F, width=12,height=12)
plotGeneNetwork(vec_kMs=list_networkPlotArgs$vec_kMs, 
                mat_adj=list_networkPlotArgs$mat_adj)
dev.off()