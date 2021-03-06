frag.relations.by.layers <- function(graph, layer.attr){
  if(! is.igraph(graph)) stop("Not an igraph object")
  if(is.null(vertex_attr(graph, layer.attr))) stop("The parameter 'layer.attr' is required.")
  if( ! is.character(layer.attr))  stop("The parameter 'layer.attr' requires a character value.")
  if( ! layer.attr %in% names(vertex_attr(graph)) ){
    stop(paste("No '", layer.attr, "' vertices attribute", sep=""))
  }
  layers <- vertex_attr(graph, layer.attr)
  
  if(is.null(V(graph)$name)){
    V(graph)$name <- 1:gorder(graph)
  }
  
  v.list <- data.frame(v = V(graph)$name, layer = layers)
  e.list <- data.frame(as_edgelist(graph))
  e.list <- merge(e.list, v.list, by.x = "X2", by.y = "v")
  e.list <- merge(e.list, v.list, by.x = "X1", by.y = "v")
  
  res <- table(e.list$layer.x, e.list$layer.y)
  diag <- diag(res)
  res <- res + matrix(res, nrow(res), byrow=TRUE)
  diag(res) <- diag
  res[upper.tri(res)] <- NA
  res
}
