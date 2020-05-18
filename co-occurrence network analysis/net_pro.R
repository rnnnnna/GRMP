net_pro<-function(igraph,outdir){
  # network property
  # The size of the graph (number of edges)
  num.edges <- length(E(igraph)) # length(curve_multiple(igraph))
  num.edges
  # Order (number of vertices) of a graph
  num.vertices <- length(V(igraph))# length(diversity(igraph, weights = NULL, vids = 	V(igraph)))
  num.vertices
  # (connectance) 
  connectance <- edge_density(igraph,loops=FALSE)
  connectance
  # (Average degree)
  average.degree <- mean(igraph::degree(igraph))
  average.degree
  # (Average path length)
  average.path.length <- average.path.length(igraph)  # mean_distance calculates the average path length in a graph
  average.path.length
  # Diameter
  diameter <- diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NA)
  diameter
  # edge connectivity / group adhesion
  edge.connectivity <- edge_connectivity(igraph)
  edge.connectivity
  #(Clustering coefficient)
  clustering.coefficient <- transitivity(igraph) 
  clustering.coefficient
  no.clusters <- no.clusters(igraph)
  no.clusters
  # (Degree centralization)
  centralization.degree <- centralization.degree(igraph)$centralization
  centralization.degree
  # (Betweenness centralization)
  centralization.betweenness <- centralization.betweenness(igraph)$centralization 
  centralization.betweenness
  # (Closeness centralization)
  centralization.closeness <- centralization.closeness(igraph)$centralization
  centralization.closeness
  
  num.pos.edges<-sum(E(igraph)$weight>0)# number of postive correlation
  num.neg.edges<-sum(E(igraph)$weight<0)# number of negative correlation
  
  igraph.network.pro <- rbind(num.edges,num.pos.edges,num.neg.edges,num.vertices,connectance,average.degree,average.path.length,diameter,edge.connectivity,clustering.coefficient,no.clusters,centralization.degree,centralization.betweenness,centralization.closeness)
  rownames(igraph.network.pro)<-c("num.edges","num.pos.edges","num.neg.edges","num.vertices","connectance","average.degree","average.path.length","diameter","edge.connectivity","clustering.coefficient","no.clusters","centralization.degree","centralization.betweenness","centralization.closeness")
  colnames(igraph.network.pro)<- "value"
  igraph.network.pro
}