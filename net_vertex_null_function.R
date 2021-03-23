# null models for networks
network <- sample_gnp(10, 2/10)
# randomizes node identity but keeps nature/number of edges
# may not be the null model we want?
metric.obs <- diameter(network)
network.adj <- as_adjacency_matrix(network)
network.adj.d <- as.dist(network.adj)
network.adj.m <- as.matrix(network.adj)
network.dim <- dim(network.adj.m)[1]
network.label.shuf <- sample(network.dim)
network.adj.m.shuf <- network.adj.m[network.label.shuf, network.label.shuf]
network.shuf <- graph_from_adjacency_matrix(network.adj.m.shuf, mode="undirected")
diameter(network)
diameter(network.shuf)
metric.rnd <- diameter(network.shuf)
g <- sample_gnp(10, 2/10)
h <- graph_from_adjacency_matrix(as_adjacency_matrix(g), mode="undirected")
diameter(g)
diameter(h)
i.d <- as.dist(as_adjacency_matrix(g))
i <- graph_from_adjacency_matrix(i.d, mode="undirected")
diameter(i)
net_vertex_null <- function (network, runs = 999, model="unrestricted")
{
  # actually randomize edges
  network.adj <- as_adjacency_matrix(network)
  network.adj.d <- as.dist(network.adj)
  network.adj.m <- as.matrix(network.adj.d)
  
  ###extract vertex names###
  if (!is.null(network.adj@Dimnames[[1]])){
    vertex.names <- network.adj@Dimnames[[1]]
  }
  else{
    vertex.names <- c(1:dim(network.adj)[1])
  }
  
  
  ###list of metrics###
  metrics.list <- c("authority_score","closeness","coreness","degree","estimated_betweenness","average_nearest_neighbor_degree","strength","betweenness_centrality","closeness_centrality","degree_centrality","eigen_centrality")
  
  #authority_score
  authority_score.obs <- authority_score(network)$vector 
  authority_score.rnd <- matrix(0,nrow = length(vertex.names))
  
  
  #closeness
  closeness.obs <- closeness(network) 
  closeness.rnd <- matrix(0,nrow = length(vertex.names))
  
  
  #coreness
  coreness.obs <- coreness(network)
  coreness.rnd <- matrix(0,nrow = length(vertex.names))
  
  #degree
  degree.obs <- degree(network) 
  degree.rnd <- matrix(0,nrow = length(vertex.names))
  
  
  
  #estimated_betweenness
  estimated_betweenness.obs <- estimate_betweenness(network,cutoff = 0) 
  estimated_betweenness.rnd <- matrix(0,nrow = length(vertex.names))
  
  
  #average_nearest_neighbor_degree
  average_nearest_neighbor_degree.obs <- knn(network)$knn 
  average_nearest_neighbor_degree.rnd <- matrix(0,nrow = length(vertex.names))
  
  
  #strength
  strength.obs <- strength(network) 
  strength.rnd <- matrix(0,nrow = length(vertex.names))
  
  
  #betweenness_centrality
  betweenness_centrality.obs <- centr_betw(network)[[1]]
  betweenness_centrality.rnd <- matrix(0,nrow = length(vertex.names))
  names(betweenness_centrality.obs) <- vertex.names
  
  
  #closeness_centrality
  closeness_centrality.obs <- centr_clo(network)[[1]]
  closeness_centrality.rnd <- matrix(0,nrow = length(vertex.names))
  names(closeness_centrality.obs) <- vertex.names
  
  
  #degree_centrality
  degree_centrality.obs <- centr_degree(network)[[1]]
  degree_centrality.rnd <- matrix(0,nrow = length(vertex.names))
  names(degree_centrality.obs) <- vertex.names
  
  
  #eigen_centrality
  eigen_centrality.obs <- centr_eigen(network)[[1]]
  eigen_centrality.rnd <- matrix(0,nrow = length(vertex.names))
  names(eigen_centrality.obs) <- vertex.names
  
  for (i in 1:runs) 
  {
    if (model=="unrestricted"){
      
      network.adj.m.shuf <- network.adj.m
      network.adj.m.shuf[lower.tri(network.adj.m.shuf)]  <- sample(network.adj.d)
      network.adj.d.shuf <- as.dist(network.adj.m.shuf)
      network.shuf <- graph_from_adjacency_matrix(network.adj.d.shuf, mode="undirected")
      vertex_attr(network.shuf)$name <- vertex_attr(network)$name
      
    }
    
    else if (model=="restricted"){
      
      network.dim <- dim(network.adj.m)[1]
      network.label.shuf <- sample(network.dim)
      network.adj.m.shuf <- network.adj.m[network.label.shuf, network.label.shuf]
      network.shuf <- graph_from_adjacency_matrix(network.adj.m.shuf, mode="undirected")
      vertex_attr(network.shuf)$name <- vertex_attr(network)$name
    }
    
    
    
    #authority_score
    authority_score.rnd <- cbind(authority_score.rnd,authority_score(network.shuf)$vector) 
    
    #closeness
    closeness.rnd <- cbind(closeness.rnd,closeness(network.shuf)) 
    
    #coreness
    coreness.rnd <- cbind(coreness.rnd,coreness(network.shuf))
    
    #degree
    degree.rnd <- cbind(degree.rnd,degree(network.shuf)) 
    
    #estimated_betweenness
    estimated_betweenness.rnd <- cbind(estimated_betweenness.rnd,estimate_betweenness(network.shuf,cutoff = 0)) 
    
    #average_nearest_neighbor_degree
    average_nearest_neighbor_degree.rnd <- cbind(average_nearest_neighbor_degree.rnd,knn(network.shuf)$knn) 
    
    #strength
    strength.rnd <- cbind(strength.rnd,strength(network.shuf)) 
    
    #betweenness_centrality
    rnd.cb <- centr_betw(network.shuf)[[1]]
    names(rnd.cb) <- vertex.names
    betweenness_centrality.rnd <- cbind(betweenness_centrality.rnd,rnd.cb)
    
    #closeness_centrality
    rnd.cc <- centr_clo(network.shuf)[[1]]
    names(rnd.cc) <- vertex.names
    closeness_centrality.rnd <- cbind(closeness_centrality.rnd,rnd.cc)
    
    #degree_centrality
    rnd.cd <- centr_degree(network.shuf)[[1]]
    names(rnd.cd) <- vertex.names
    degree_centrality.rnd <- cbind(degree_centrality.rnd,rnd.cd)
    
    #eigen_centrality
    rnd.ce <- centr_eigen(network.shuf)[[1]]
    
    names(rnd.ce) <- vertex.names
    
    eigen_centrality.rnd <- cbind(eigen_centrality.rnd,rnd.ce)
    
  }
  
  results.list <- vector("list", 11)
  #iterate over each metric and calculate statistics
  for (n in c(1:length(metrics.list))) {
    metric <- metrics.list[n]
    results <- data.frame()
    
    metric.dat <- as.data.frame(get(paste0(metric,".rnd")))
    
    rownames(metric.dat) <- vertex.names
    
    
    
    for (vertex in vertex.names){
      
      metric.obs <- get(paste0(metric,".obs"))[vertex]
      metric.rnd <- as.matrix(metric.dat[vertex,])
      
      
      metric.rand.mean <- mean(metric.rnd)
      
      metric.rand.sd <- sd(metric.rnd)
      
      metric.obs.z <- (metric.obs - metric.rand.mean) / metric.rand.sd
      
      metric.obs.rank <- rank(c(metric.obs,metric.rnd))[1]
      metric.obs.p <- metric.obs.rank / (runs + 1)
      
      stats <- c(metric.obs, metric.rand.mean, metric.rand.sd, metric.obs.z, metric.obs.rank, metric.obs.p)
      results <- rbind(results,stats)
      
      
    }
    
    #put the vertex names
    rownames(results) <- vertex.names
    colnames(results) <- c("Observed metric", "Null model mean", "Null model standard deviation", "Observed metric z score", "Observed metric rank", "Observed metric p-value")
    
    results.list[n] <- list(results)
    
    
    
    
    
  }
  
  names(results.list) <- metrics.list
  
  return(results.list)
  
  
}
# example
test1 <- net_vertex_null(g)
test2 <- net_vertex_null(g,model = "restricted")
test3 <- net_vertex_null(ig.sparcc)
test4 <- net_vertex_null(ig.sparcc,model = "restricted")
vertex_connectivity(g)