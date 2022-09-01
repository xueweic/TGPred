





#' Estimate selection probability using MSENet function solving by APGD
#'
#' @param X expressional levels of n_genes target genes (TGs)
#' @param y expressional levels of a transcription factor (TF)
#' @param Adj the adjacency matrix of network structure.
#' @param alphas the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
#' @param n_lambda the number of lambdas
#' @param ratio the ratio of smallest lambda. default: 0.01
#' @param B the number of half-sample resamplings used to calculate selection probabilities of genes. default: 500
#' @param gamma initial value of gamma in APGD. default: 1000
#' @param niter the maximum number of APGD to solve Net regression. default: 2000
#' @param crit_beta converge criterion of change of beta. default: 1e-4
#' @param crit_obj converge criterion of change of objective function. default: 1e-8
#' @param timer decide if exist the output report. default: FALSE
#'
#' @return SP: n_genes length vector of selection probability.
#' @export
#'
#' @examples
MSENet_SP <- function(X, y, Adj, alphas, n_lambda, ratio=1e-2, B=500,  gamma=1000, niter=2000,
                        crit_beta=1e-4, crit_obj=1e-8, timer=TRUE){
  X.ori <- data.matrix(X)
  X <- scale(X.ori)
  y.ori <- data.matrix(y)
  y <- scale(y.ori)
  # ---- check X y
  if (nrow(X) != nrow(y)){
    stop("Error: Please check the sample size of X and y. They should be the same!")
  } else if (ncol(y) != 1){
    stop("Error: Please check the dimension of y. It should be n*1 vector!")
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
  # ---- check alphas
  if (sum(alphas <= 0) + sum(alphas > 1) != 0){
    stop("Error: The range of alpha shoud be in (0,1]!")
  }
  # if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
  library(glmnet)
  ####### Setting:
  n <- nrow(X)
  p <- ncol(X)
  # - Laplacian matrix
  res <- CalculateLaplacian(Adj)
  L <- res$L
  L_norm <- res$L_norm
  # Calculate S: ridge regression:
  lambda_seq <- 10^seq(1, -2, by = -.01)
  ridge_cv <- cv.glmnet(X, y, alpha = 0, lambda = lambda_seq)
  # Best lambda value
  best_lambda <- ridge_cv$lambda.min
  best_ridge <- glmnet(X, y, alpha = 0, lambda  = best_lambda)
  beta_hat <- best_ridge$beta@x
  S <- diag(sign(beta_hat))
  SLS <- L_norm * diag(S)
  SLS <- t(SLS) * diag(S)

  ## Start MSENet by APGD method
  print("Start calculating selection probability using MSENet by APGD method:")
  n_alpha <- length(alphas)
  SP.LambdaAlpha <- matrix(NA, nrow = p, ncol = n_lambda*n_alpha)
  flag <- 0
  start.time <- proc.time()
  for (i.alpha in 1:n_alpha){
    alpha <- alphas[i.alpha]

    lambda_set <- Lambda_grid(X.ori, y.ori, n_lambda, alpha, loss_func = "MSE", ratio)
    for (i.lambda in 1:length(lambda_set)){
      lambda <- lambda_set[i.lambda]
      flag <- flag + 1

      print(paste0("Case ", flag, " : lambda ", lambda, ", alpha ", alpha))
      beta_hat.LambdaAlpha <- matrix(NA, nrow = p, ncol = B)
      for (i.b in 1:B){
        pos <- sample(1:n, n/2)
        y.sub <- y[pos,]
        X.sub <- X[pos,]
        invisible(capture.output(beta_hat_APGD <- MSENet_Beta(X.sub, y.sub, Adj, lambda, alpha, method="APGD",
                                  gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=TRUE)))
        beta_hat.LambdaAlpha[,i.b] <- beta_hat_APGD
      }
      Proportion.LambdaAlpha <- (beta_hat.LambdaAlpha != 0) * 1
      SP.LambdaAlpha[,flag] <- rowMeans(Proportion.LambdaAlpha)
    }
  }
  SP.Net <- apply(SP.LambdaAlpha, 1, max)
  print("Done with MSENet by APGD method!")
  end.time <- proc.time()
  if (timer == TRUE){
    diff <- end.time - start.time
    cost_h <- round(diff[3]/3600, 4)
    cost_min <- round(diff[3]/60, 4)
    cost_s <- round(diff[3], 4)
    print(paste("The calculate time:", cost_h, "h", cost_min, "min", cost_s, "s."))
  }
  return(SP.Net)
}





