
#' Simulate a covariance matrix from the specific graph (Adj) based on a Gaussian graphical model.
#'
#' @param Adj the adjacency matrix of network structure
#' @param a1 parameters for constructing domain [a1, a2] union [b1, b2]
# 					default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7
#' @param a2 parameters for constructing domain [a1, a2] union [b1, b2]
# 					default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7
#' @param b1 parameters for constructing domain [a1, a2] union [b1, b2]
# 					default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7
#' @param b2 parameters for constructing domain [a1, a2] union [b1, b2]
# 					default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7
#'
#' @return the covariance matrix of genes according to network structure.
#' @export
#'
#' @examples


GraphicalModel <- function(Adj, a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7) {
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  library(MASS)
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

  ## Calculate Sigma
  n <- nrow(Adj) * ncol(Adj)
  rho <- runif(n, 0, a2 - a1 + b2 - b1)
  rho.1 <- (a1 + rho) * (rho < (a2 - a1)) + (b1 + rho - (a2 - a1)) * (rho >= (a2 - a1))
  Adj <- Adj * matrix(rho.1, nrow(Adj))
  # - Rescale matrix
  rSum <- rowSums(abs(Adj))
  Adj.rescale <- Adj / (5 * sapply(1:ncol(Adj), function(t) {
    return(rSum)
  }))
  # - Ensure symmetry
  A <- 0.5 * (Adj.rescale + t(Adj.rescale))
  diag(A) <- 1
  # - Calculate Sigma
  A.inv <- ginv(A)
  A.diag <- as.matrix(diag(A.inv))
  Sigma <- A.inv / sqrt(A.diag %*% t(A.diag))
  return(Sigma)
}
