

#' Simulate a grid set of lambdas for a given alpha in penalized regression.
#'
#' @param X expression levels of n_genes target genes (TGs).
#' @param y expression levels of a transcription factor (TF).
#' @param n_lambda the number of lambdas.
#' @param alpha the proportion of l1 norm affects (the numerical values of nonzero coefficients). It's in range (0,1].
#' @param loss_func either "Huber" or "MSE".
# If loss_func = "Huber", the loss function in penalized regression model is Huber function.
# If loss_func = "MSE", the loss function in penalized regression model is mean squared errors.
#'
#' @return n_lambda length vector of lambdas according to the alpha you provided.
#' @export
#'
#' @examples
Lambda_grid <- function(X, y, n_lambda, alpha, loss_func) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  ## Check X & y
  if (nrow(X) != nrow(y)) {
    stop("Error: Please check the sample size of y and X. They should be the same!")
  } else if (ncol(y) != 1) {
    stop("Errro: Please check the number of y. It should be 1!")
  }
  ##  Check alpha range
  if (alpha <= 0 || alpha > 1) {
    stop("Error: The range of alpha should be in (0,1]!")
  }
  X <- scale(X)
  y <- scale(y)
  ratio <- 0.01
  ## Check loss function
  # Use Huber function
  if (loss_func == "Huber") {
    M <- IQR(t(X) %*% y) / 1.345
    H_grad <- Huber_grid(y, M)
    lambda_max <- max(t(H_grad) %*% X) / alpha
  }
  # Use mean squared errors function
  else if (loss_func == "MSE") {
    lambda_max <- max(t(y) %*% X) / alpha
  }
  # Not defined
  else {
    stop("Error: Please set the loss funtion from 'Huber' or 'MSE'!")
  }
  lambdas <- exp(log(10) * seq(log10(lambda_max * ratio), log10(lambda_max), length.out = n_lambda))
  return(lambdas)
}
