library(ggplot2)
requireNamespace("igraph")
requireNamespace("reshape")


filtrationCurves <- function(nets, nThres, desc, p1){
  # filtration: native edge weights
  thres <- quantile(unlist(nets), seq(0, 1, length.out = nThres))
  thres <- rev(thres) # growing from the largest weights

  # graph descriptor on each subgraph
  fc <- matrix(0, nrow = length(nets), ncol = nThres)
  if (desc == 'X' | desc == 'Y'){ # node label count
    for (i in 1:nThres){
      print(i)
      # node label counts with filtered edges
      nlc <- lapply(nets, function(net){
        subnet <- net > thres[i]
        labels <- union(which(rowSums(subnet) > 0), which(colSums(subnet) > 0))
        count <- length(labels)
        return(count)
      })
      # node label histogram
      fc[, i] <- unlist(nlc)
    }
  } else if (grepl('lcc', desc)){ # largest connected component
    for (i in 1:nThres){
      print(i)
      fc_i <- lapply(nets, function(net){
        net <- ifelse(net >= thres[i], 1, 0)
        g <- igraph::graph_from_adjacency_matrix(net, mode = 'undirected')
        lcc <- igraph::components(g)
        if (desc == 'lccs'){
          return(max(lcc$csize))
        } else if (desc == 'lccd'){
          lccIdx <- which.max(lcc$csize)
          lccIds <- igraph::V(g)[lcc$membership == lccIdx]
          gsub <- igraph::induced_subgraph(g, lccIds)
          return(mean(igraph::degree(gsub, V(gsub), normalized = F)))
        }
      })
      fc[, i] <- unlist(fc_i)
    }
  }
  return(list(thres = thres, fc = fc))
}


plotCurves <- function(fc, desc){
  thres <- seq(0, 1, length.out=length(fc$thres))
  fc <- fc$fc
  if (desc == 'X'){
    c1 <- 'lightblue'
    c2 <- 'blue'
    ylabel <- 'Node label count (X)'
    } else if (desc == 'Y'){
    c1 <- 'lightsalmon'
    c2 <- 'red'
    ylabel <- 'Node label count (Y)'
    } else if (label == 'lcc'){
      c1 <- 'lightblue'
      c2 <- 'blue'
      ylabel <- 'Mean degree'
    }
  df <- data.frame(weight=thres, t(fc))
  df <- reshape::melt(df, id='weight', variable_name='ISN')
  
  ggplot() +
    geom_line(data=df, aes(x=weight, y=value, group=ISN), 
              color=c1, size = 0.5, alpha=0.5) +
    geom_line(aes(x=thres, y=colMeans(fc)), 
              color=c2, size = 0.5, alpha=1) +
    scale_x_continuous(labels = scales::percent, breaks = thres) +
    theme_linedraw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15)) +
    ggtitle(paste0('Filtration curves of ISNs (largest connected component)')) +
    xlab('Edge weight threshold') + ylab(ylabel)
}


plotCurvesClust <- function(fc, clust, desc){
  thres <- seq(0, 1, length.out=length(fc$thres))
  fc <- fc$fc
  if (desc == 'X'){
    ylabel <- 'Node label count (X)'
  } else if (desc == 'Y'){
    ylabel <- 'Node label count (Y)'
  } else if (desc == 'lcc'){
    ylabel <- 'Mean degree'
  }
  
  clustLabel <- unique(clust)
  clusters <- lapply(clustLabel, function(c){which(clust == c, arr.ind = T)})
  clustMean <- c()
  clustSD <- c()
  for (cluster in clusters){
    clustMean <- c(clustMean, colMeans(fc[cluster,]))
    clustSD <- c(clustSD, apply(fc[cluster,], 2, sd))
  }
  df <- data.frame(cluster=rep(clustLabel, each=length(thres)), 
                  weight=rep(thres, length(clusters)),
                  mean=clustMean, sd=clustSD)
  
  pd <- position_dodge(width = 0.02)
  ggplot(data=df, aes(x=weight, y=mean, group=cluster, color=cluster)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.05, position=pd) +
    geom_line() +
    scale_x_continuous(labels = scales::percent_format(accuracy = 5L), breaks = thres) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15),
          legend.title = element_text(size=20),
          legend.text = element_text(size=15),
          legend.position = "bottom") +
    ggtitle(paste0('Filtration curves of ISNs in every cluster')) +
    xlab('Edge weight threshold') + ylab(ylabel)
}




