


#' Simulate y and X from a given network structure.
#'
#' @param N_samples the number of sample size.
#' @param N_genes the number of traget genes.
#' @param Adj the adjacency matrix of network structure.
#' Adjacency matrix must be a N_genes*N_genes dimensional symmetric matrix,
#' the elements equal 1 indicates two genes are connected.
#' If you consider Barabasi-Albert Network or Hierarchical Network in the article,
#' you can directly use "ConstructNetwork" function to get the adjacency matrix.
#' @param Sigma the covariance matrix of target genes according to network structure.
#' You can directly use "GraphicalModel" function to get the covariance matrix.
#' @param method "HN": by Hierarchical Network, "BAN": by Barabasi-Albert Network or "DIY": by user designed
#' @param beta0 numeric value of effect size in simulation settings.
#' 		   default: NULL; if method is "HN" or "BAN", input a nunerical value.
#' @param beta_true numeric matrix with the dimension of N_genes * 1 in simulation settings.
#' 	   default: NULL; if method is "DIY", input a nunerical matrix (N_genes * 1).
#'
#' @return a list of y: expression levels of a transcription factor (TF);
#' X: expression levels of n_genes target genes (TGs);
#' beta: true regulated effect beta for N_genes TGs.
#' @export
#'
#' @examples

SimulationData <- function(N_samples, N_genes, Adj, Sigma, method, beta0 = NULL, beta_true = NULL) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) install.packages("mvtnorm")
  library(mvtnorm)
  ## Check method:
  if (method == "DIY") {
    if (!is.null(beta_true)) {
      if (length(beta_true) != N_genes) {
        stop("Error: In 'DIY' method, beta_true must be a numerical vector with the dimension of N_genes!")
      }
    } else {
      stop("Error: In 'DIY' method, please give a numerical vector with dimensiona of N_genes to beta_true!")
    }
    beta_true <- as.matrix(beta_true)
  } else if (method == "HN" || method == "BAN") {
    if (!is.null(beta0)) {
      if (length(beta0) != 1) {
        stop("Error: In 'HN' and 'BAN' method, beta0 must be a numerical value!")
      }
    } else {
      stop("Error: In 'HN' and 'BAN' methods, please give a numerical value to beta0 (ex. beta0=0.1)!")
    }
  } else {
    stop("Error: Please check method, must be one of HN, BAN, or DIY!")
  }
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
  ## Check Sigma
  if (nrow(Sigma) != ncol(Sigma)) {
    stop("Error: The covariance matrix must be squared matrix! Please try to obtain
             the covariance matrix Sigma from adjacency matrix using 'GraphicalModel' function!")
  } else {
    if (!all.equal(Sigma, t(Sigma))) {
      stop("Error: The covariance matrix must be symmetric matrix! Please try to obtain
             the covariance matrix Sigma from adjacency matrix using 'GraphicalModel' function!")
    } else if (unique(diag(Sigma)) != 1) {
      stop("Error: The diagnal elements of covariance matrix must be all 1! Please try to obtain
             the covariance matrix Sigma from adjacency matrix using 'GraphicalModel' function!")
    }
  }


  # - step 1: generate y
  y <- as.matrix(rnorm(N_samples, mean = 0, sd = 1))
  # - step 2: calculate degree
  n_degree <- rowSums(Adj)
  # - step 3: generate X
  beta <- rep(0, N_genes)
  if (method == "HN") {
    beta[1] <- beta0
    geneAffact.1 <- c(seq(2, 12), seq(24, 34))
    beta[geneAffact.1] <- sqrt(n_degree[geneAffact.1]) * beta0 / 3
    geneAffact.2 <- c(seq(13, 23), seq(35, 45))
    beta[geneAffact.2] <- -sqrt(n_degree[geneAffact.2]) * beta0 / 3
  } else if (method == "BAN") {
    pos <- c(1:10, 21:30)
    beta[pos] <- sqrt(n_degree[pos]) * beta0
    pos <- c(11:20, 31:40)
    beta[pos] <- -sqrt(n_degree[pos]) * beta0
  } else if (method == "DIY") {
    beta <- beta_true
  }
  beta <- as.matrix(beta)
  Z <- rmvnorm(N_samples, mean = rep(0, N_genes), sigma = Sigma)
  error <- rmvnorm(N_samples, mean = rep(0, N_genes), sigma = diag(rep(1, N_genes)))
  X <- y %*% t(beta) + Z + error
  return(list("y" = y, "X" = X, "beta" = beta))
}
