





#' Estimate selection probability using Lasso function solving by APGD
#'
#' @param X expressional levels of n_genes target genes (TGs)
#' @param y expressional levels of a transcription factor (TF)
#' @param n_lambda the number of lambdas
#' @param B the number of half-sample resamplings used to calculate selection probabilities of genes. default: 500
#' @param gamma initial value of gamma in APGD. default: 1000
#' @param niter the maximum number of APGD to solve Lasso regression. default: 2000
#' @param crit_beta converge criterion of change of beta. default: 1e-4
#' @param crit_obj converge criterion of change of objective function. default: 1e-8
#' @param timer decide if exist the output report. default: FALSE
#'
#' @return SP: n_genes length vector of selection probability.
#' @export
#'
#' @examples
Lasso_SP <- function(X, y, n_lambda, B=500,  gamma=1000, niter=2000,
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
  # if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
  library(glmnet)
  ####### Setting:
  n <- nrow(X)
  p <- ncol(X)


  ## Start Huber-Net by APGD method
  print("Start calculating selection probability using Lasso by APGD method:")
  SP.LambdaAlpha <- matrix(NA, nrow = p, ncol = n_lambda)
  flag <- 0
  start.time <- proc.time()

  lambda_set <- Lambda_grid(X.ori, y.ori, n_lambda, alpha=1, loss_func = "MSE")
  for (i.lambda in 1:length(lambda_set)){
    lambda <- lambda_set[i.lambda]
    flag <- flag + 1

    print(paste0("Case ", flag, " : lambda ", lambda))
    beta_hat.LambdaAlpha <- matrix(NA, nrow = p, ncol = B)
    for (i.b in 1:B){
      pos <- sample(1:n, n/2)
      y.sub <- y[pos,]
      X.sub <- X[pos,]
      beta_hat_APGD <- Lasso_Beta(X, y, lambda, method="APGD",
                                  gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=TRUE)
      beta_hat.LambdaAlpha[,i.b] <- beta_hat_APGD
    }
    Proportion.LambdaAlpha <- (beta_hat.LambdaAlpha != 0) * 1
    SP.LambdaAlpha[,flag] <- rowMeans(Proportion.LambdaAlpha)
  }

  SP.Lasso <- apply(SP.LambdaAlpha, 1, max)
  print("Done with Lasso by APGD method!")
  end.time <- proc.time()
  if (timer == TRUE){
    diff <- end.time - start.time
    cost_h <- round(diff[3]/3600, 4)
    cost_min <- round(diff[3]/60, 4)
    cost_s <- round(diff[3], 4)
    print(paste("The calculate time:", cost_h, "h", cost_min, "min", cost_s, "s."))
  }
  return(SP.Lasso)
}
