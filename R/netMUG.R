library(parallel)
library(SmCCNet)
library(devtools)
library(dynamicTreeCut)


selectFeatures <- function(X, Y, Z, l1, l2, s1, s2, 
                              noTrait=FALSE, subsamplingNum=500, CCcoef=NULL,
                              cutHeight=0.999, plotTree=FALSE){
  
  if(length(dim(X)) != 2){
    stop("The current version of netMUG requires X to be a 2D matrix or data frame.")
  }
  if(length(dim(Y)) != 2){
    stop("The current version of netMUG requires Y to be a 2D matrix or data frame.")
  }
  if(dim(X)[1] != dim(Y)[1]){
    stop("The first dimension of X and Y doesn't match.")
  }
  
  n <- dim(X)[1]; p1 <- dim(X)[2]; p2 <- dim(Y)[2]
  
  if(is.null(colnames(X))){
    colnames(X) <- paste0("X_", seq_len(p1))
  }
  if(is.null(colnames(Y))){
    colnames(Y) <- paste0("Y_", seq_len(p2))
  }
  
  # Step 1: Calculate the canonical correlation weights based on SmCCA
  Ws <- getRobustPseudoWeights(X, Y, Trait = Z, Lambda1 = l1, Lambda2 = l2, 
                               s1 = s1, s2 = s2, NoTrait = noTrait, FilterByTrait = FALSE,
                               SubsamplingNum = subsamplingNum, CCcoef=CCcoef, trace = FALSE)
  Ws_nonzero <- which(rowSums(abs(Ws)) != 0)
  
  # Step 2: Compute the similarity matrix based on one or more canonical correlation weight vectors
  featureLabel <- c(colnames(X), colnames(Y))
  sim <- getAbar(Ws[Ws_nonzero,], 
                 P1 = which(Ws_nonzero > p1)[1]-1, 
                 FeatureLabel = featureLabel[Ws_nonzero])
  
  # Step 3: Extract multi-view modules based on the similarity matrix
  modules <- getMultiOmicsModules(sim, 
                                 P1 = tail(which(Ws_nonzero <= p1), 1), 
                                 CutHeight = cutHeight, PlotTree = plotTree)
  #modules <- lapply(modules, function(module) Ws_nonzero[module])
  featureIdx <- Ws_nonzero[unlist(modules)]
  featureX <- featureIdx[featureIdx <= p1]
  featureY <- featureIdx[featureIdx > p1] - p1
  return(list(featureX = featureX, featureY = featureY,
              Ws = Ws, sim = sim, modules = modules))
}


plotModules <- function(X, Y, Ws, sim, modules, corr=NULL, moduleIdx=1,
                        edgeCut=0, addCorrSign=TRUE, saveFile=NULL, 
                        showType1Label=TRUE, showType2Label=TRUE, 
                        plotTitle="", netLayout="lgl",
                        showNodes=TRUE, vertexLabelCex=1, vertexSize=1){
  
  Ws_nonzero <- which(rowSums(abs(Ws)) != 0)
  if(is.null(corr)){
    corr <- cor(cbind(X, Y)[,Ws_nonzero])
  }
  featureLabel <- c(colnames(X), colnames(Y))[Ws_nonzero]
  p1 <- which(Ws_nonzero > dim(X)[2])[1]-1
  
  plotMultiOmicsNetwork(sim, corr, modules, ModuleIdx = moduleIdx,
                        P1 = p1, AddCorrSign = addCorrSign, EdgeCut = edgeCut,
                        FeatureLabel = featureLabel, SaveFile = saveFile,
                        ShowType1Label = showType1Label, ShowType2Label = showType2Label,
                        NetLayout = netLayout, ShowNodes = showNodes, 
                        VertexSize = vertexSize, VertexLabelCex = vertexLabelCex)
}


buildInfISNs <- function(V, Z, nCores=1){
  n <- dim(V)[1]; p <- dim(V)[2]
  gCorV <- cor(V)
  gCorZ <- cor(V, Z)
  gnet <- gCorV + matrix(gCorZ, nrow=p, ncol=p, byrow=TRUE) + matrix(gCorZ, nrow=p, ncol=p, byrow=FALSE)
  cl <- makeCluster(nCores, type = "PSOCK")
  clusterExport(cl, c("V", "Z", "gnet", "n", "p"))
  ISNs <- parLapply(cl, seq_len(n), function(i){
    corV <- cor(V[-i,])
    corZ <- cor(V[-i,], Z[-i])
    net <- corV + matrix(corZ, nrow=p, ncol=p, byrow=TRUE) + matrix(corZ, nrow=p, ncol=p, byrow=FALSE)
    ISN <- abs(gnet - net)
    diag(ISN) <- 0
    return(ISN)
  })
  stopCluster(cl)
  return(ISNs)
}


computeDist <- function(ISNs){
  p <- dim(ISNs[[1]])[1]
  ind <- expand.grid(1:p, 1:p)
  ind <- as.matrix(ind[which(ind[,1] > ind[,2]), ])
  edges <- do.call(rbind, lapply(ISNs, function(net){as.numeric(net[ind])}))
  distEuc <- dist(scale(edges))
  distEuc <- as.matrix(distEuc)
  return(distEuc)
}


netMUG <- function(X, Y, Z, l1, l2, s1, s2, noTrait=FALSE, subsamplingNum=500, 
                   CCcoef=NULL, cutHeight=0.999, plotTree=FALSE, nCores=1,
                   minClusterSize=1, deepSplit=0){
  
  # Step 1: select multi-view features informed by an extraneous variable
  smccnet <- selectFeatures(X, Y, Z, l1, l2, s1, s2, noTrait=noTrait, 
                           subsamplingNum=subsamplingNum, CCcoef=CCcoef,
                           cutHeight=cutHeight, plotTree=plotTree)
  Xsub <- X[, smccnet$featureX]
  Ysub <- Y[, smccnet$featureY]
  
  # Step 2: build ISNs from the selected features
  V <- cbind(Xsub, Ysub)
  ISNs <- buildInfISNs(V, Z, nCores = nCores)
  
  # Step 3: compute distances between ISNs
  dis <- computeDist(ISNs)
  
  # Step 4: Ward's hierarchical clustering with Dynamic Tree Cut
  dendro <- hclust(as.dist(dis_m), method = "ward.D2")
  clust <- cutreeDynamic(dendro, minClusterSize = minClusterSize, distM = dis, 
                        deepSplit = deepSplit)
  clust <- as.factor(clust)
  
  return(list(Xsub = Xsub, Ysub = Ysub, ISNs = ISNs, clust = clust))
}





