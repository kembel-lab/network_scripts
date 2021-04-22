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
net_null <- function (network, runs = 999)
{
  # actually randomize edges
  network.adj <- as_adjacency_matrix(network)
  network.adj.d <- as.dist(network.adj)
  network.adj.m <- as.matrix(network.adj.d)
  
  ###list of metrics###
  metrics.list <- c("diameter","motifs","edge_connectivity","vertex_connectivity","girth","radius","centralized_betweenness","centralized_closeness","centralized_degree","centralized_eigen_centrality")
  # diameter
  diameter.obs <- diameter(network)
  diameter.rnd <- vector(length=runs)
  # number of motifs
  motifs.obs <- count_motifs(network)
  motifs.rnd <- vector(length=runs)
  # edge_connectivity
  edge_connectivity.obs <- edge_connectivity(network)
  edge_connectivity.rnd <- vector(length=runs)
  # vertex_connectivity
  vertex_connectivity.obs <- vertex_connectivity(network)
  vertex_connectivity.rnd <- vector(length=runs)
  # girth
  girth.obs <- girth(network)[[1]]
  girth.rnd <- vector(length=runs)
  # radius
  radius.obs <- radius(network)
  radius.rnd <- vector(length=runs)
  # centralized_betweenness
  centralized_betweenness.obs <- centr_betw(network)$centralization
  centralized_betweenness.rnd <- vector(length=runs)

  # centralized_closeness
  centralized_closeness.obs <- centr_clo(network)$centralization
  centralized_closeness.rnd <- vector(length=runs)

  # centralized_degree
  centralized_degree.obs <- centr_degree(network)$centralization
  centralized_degree.rnd <- vector(length=runs)

  # centralized_eigen_centrality
  centralized_eigen_centrality.obs <- centr_eigen(network)$centralization
  centralized_eigen_centrality.rnd <- vector(length=runs)

  
    for (i in 1:runs) 
  {
    network.adj.m.shuf <- network.adj.m
    network.adj.m.shuf[lower.tri(network.adj.m.shuf)]  <- sample(network.adj.d)
    network.adj.d.shuf <- as.dist(network.adj.m.shuf)
    network.shuf <- graph_from_adjacency_matrix(network.adj.d.shuf, mode="undirected")
    
    # to change metric, change this line
    #metric.rnd[i] <- diameter(network.shuf)
  
    # diameter
    diameter.rnd[i] <- diameter(network.shuf)
    
    # number of motifs
    motifs.rnd[i] <- count_motifs(network.shuf)
    
    # edge_connectivity
    edge_connectivity.rnd[i] <- edge_connectivity(network.shuf)
    
    # vertex_connectivity
    vertex_connectivity.rnd[i] <- vertex_connectivity(network.shuf)
    
    # girth
    girth.rnd[i] <- girth(network.shuf)[[1]]
    
    # radius
    radius.rnd[i] <- radius(network.shuf)
    
    # centralized_betweenness
    centralized_betweenness.rnd[i] <- centr_betw(network.shuf)$centralization
    
    # centralized_closeness
    centralized_closeness.rnd[i] <- centr_clo(network.shuf)$centralization
    
    # centralized_degree
    centralized_degree.rnd[i] <- centr_degree(network.shuf)$centralization
    
    # centralized_eigen_centrality
    centralized_eigen_centrality.rnd[i] <- centr_eigen(network.shuf)$centralization
    
    
    
    
  }
  
  #iterate over each metric and calculate statistics
  results <- data.frame()
  distribution <- data.frame(matrix(nrow=runs,ncol=length(metrics.list)))
  rownames(distribution) <- c(1:runs)
  colnames(distribution) <- metrics.list
  
  
  for (metric in metrics.list) {
    
    
    
    metric.obs <- get(paste0(metric,".obs"))
    
    metric.rnd <- get(paste0(metric,".rnd"))
    distribution[,metric] <- metric.rnd
    metric.rand.mean <- mean(metric.rnd)
    metric.rand.sd <- sd(metric.rnd)
    metric.obs.z <- (metric.obs - metric.rand.mean) / metric.rand.sd
    metric.obs.rank <- rank(c(metric.obs,metric.rnd))[1]
    metric.obs.p <- metric.obs.rank / (runs + 1)
    
    stats <- c(metric.obs, metric.rand.mean, metric.rand.sd, metric.obs.z, metric.obs.rank, metric.obs.p)
    
    results <- rbind(results,stats)
    
    
    
  }
  rownames(results) <- metrics.list
  colnames(results) <- c("Observed metric", "Null model mean", "Null model standard deviation", "Observed metric z score", "Observed metric rank", "Observed metric p-value")
  
  return(list(results,distribution))
  
  
}
# example
test <- net_null(g)
test <- net_null(ig.sparcc)
test
vertex_connectivity(g)
