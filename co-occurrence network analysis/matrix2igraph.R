matrix2igraph<-function(matr,r.threshold,p.threshold){
  occor<-corAndPvalue(matr,method = c( "spearman"))  
  # multiple test the p values
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="BH")
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  ##R values
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0 
  
  # Construct igraph object
  igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
  igraph
  bad.vs <- V(igraph)[degree(igraph) == 0]
  igraph <- delete.vertices(igraph, bad.vs)
  igraph
}