
#' Construct network in simualtion studies
#'
#' Construct the network structure from either Hierarchical Network or Barabasi-Albert Network in simualtion studies.
#'
#' @param N_genes the number of genes
#' @param Structure "HN": Hierarchical Network or  "BAN": Barabasi-Albert Network
#'
#' @return the adjacency matrix of network structure
#' @export
#'
#' @examples Adj = ConstructNetwork(N_genes = 200, Structure = "BAN")


ConstructNetwork <- function(N_genes, Structure) {
  if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
  library(Matrix)
  library(igraph)
  if (Structure == "HN") {
    n.each <- 100
    if ((N_genes %% n.each) != 0) {
      stop("Error: In Hierarchical Network, number of genes must be the integer times 100!")
    } else {
      n.module <- N_genes / n.each
      # Adjancency matrix in each module
      Adj <- matrix(0, 100, 100)
      Adj[1, seq(2, 100, 11)] <- 1
      for (i.g in 1:9) {
        Adj.group <- matrix(0, 11, 11)
        Adj.group[1, c(2, 3, 4, 5)] <- 1
        Adj.group[2, 6] <- Adj.group[3, c(7, 8)] <- Adj.group[4, c(9, 10)] <- Adj.group[5, 11] <- 1
        n.block <- seq(11 * (i.g - 1) + 2, 11 * i.g + 1)
        Adj[n.block, n.block] <- Adj.group
      }
      Adj <- Adj + t(Adj)
      # Adjancency matrix for all network
      Adj.1 <- matrix(0, 100 * n.module, 100 * n.module)
      for (i.module in 1:n.module) {
        n.block <- seq(100 * (i.module - 1) + 1, 100 * i.module)
        Adj.1[n.block, n.block] <- Adj
      }
      Adj <- Adj.1
    }
  } else if (Structure == "BAN") {
    n.each <- 10
    if ((N_genes %% n.each) != 0) {
      stop("Error: In Barabasi-Albert Network, number of genes must be the integer times 10!")
    } else {
      n.module <- N_genes / n.each
      a <- lapply(1:n.module, function(t) {
        temp <- barabasi.game(10, directed = FALSE)
        return(data.matrix(as_adjacency_matrix(temp, type = "both")))
      })
      Adj <- data.matrix(bdiag(a))
    }
  } else {
    stop("Error: Please check 'Structure', must be 'HN' as Hierarchical Network or 'BAN' as Barabasi-Albert Network!")
  }
  return(Adj)
}
