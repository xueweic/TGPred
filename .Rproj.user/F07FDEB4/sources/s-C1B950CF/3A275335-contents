

#' Calculate adjacency matrix from an annotation file.
#'
#' @param Annotation n_genes*n_pathways dimensional matrix that indicates the annotation of genes within pathways information.
#'
#' @return the adjacency matrix of network structure.
#' @export
#'
#' @examples
CalculateAdj <- function(Annotation) {
  Anno <- data.matrix(Annotation)
  # - Adjacency matrix
  Adj <- (Anno %*% t(Anno) != 0) * 1
  diag(Adj) <- 0
  no_connection <- which(colSums(Adj) == 0)
  if (length(no_connection) != 0){
    stop(paste("Error: There are some PWGs that have no connection with other PWGs. Please remove them in your input data files (Annotation and PWG expression file). The no connected PWGs index are listed below: \n", list(no_connection)))
  }
  return(Adj)
}
