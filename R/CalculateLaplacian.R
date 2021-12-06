

#' Calculate Laplacian matrix and symmetric normalized Laplacian matrix from an adjacency matrix.
#'
#' @param Adj the adjacency matrix of network structure.
#'
#' @return a list of Laplacian matrix (L) and symmetric normalized Laplacian matrix (L_norm) for network structure.
#' @export
#'
#' @examples
CalculateLaplacian <- function(Adj) {
  ## Check adjacency matrix
  if (nrow(Adj) != ncol(Adj)) {
    stop("Error: The adjacency matrix must be squared matrix!")
  } else {
    if (!all.equal(Adj, t(Adj))) {
      stop("Error: The adjacency matrix must be symmetric matrix!")
    } else if (unique(diag(Adj)) != 0) {
      stop("Error: The diagnal elements of adjacency matrix must be all 0!")
    } else if (!all.equal(unique(as.vector(Adj)), c(0, 1))) {
      stop("Error: The elements of adjacency matrix must be 0 or 1!")
    }
  }

  # Degree matrix
  D <- diag(rowSums(Adj))
  # Laplacian matrix
  L <- D - Adj
  L_norm <- L / sqrt(diag(D))
  L_norm <- t(L_norm) / sqrt(diag(D))

  return(list("L" = L, "L_norm" = L_norm))
}
