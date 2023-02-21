


#' Estimate beta_hat using Huber Elastic Net function
#'
#' @param X expressional levels of n_genes target genes (TGs)
#' @param y expressional levels of a transcription factor (TF)
#' @param lambda0 one of parameters in Huber Elastic Net regression, which controls the number of nonzero coefficients.
#' @param alpha0 one of parameters in Huber Elastic Net regression, which controls the numerical values of nonzero coefficients.
#' @param method The current methods must be 'APGD' or 'CVX'
#' @param gamma initial value of gamma in APGD. default: 1000
#' @param niter the maximum number of APGD to solve Huber Elastic Net regression. default: 2000
#' @param crit_beta converge criterion of change of beta. default: 1e-4
#' @param crit_obj converge criterion of change of objective function. default: 1e-8
#' @param quiet decide if exist the output report. default: FALSE
#' @param if.scale decide if scale the expression levels. default: FALSE
#'
#' @return beta: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in HuberENET regression.
#' @export
#'
#' @examples
#'
HuberENET_Beta <- function(X, y, lambda0, alpha0, method="APGD",
                          gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=FALSE, if.scale=FALSE){
  X <- data.matrix(X)
  y <- data.matrix(y)
  if (if.scale == "TRUE"){
    X <- scale(X)
    y <- scale(y)
  }
  # ---- check X y
  if (nrow(X) != nrow(y)){
    stop("Error: Please check the sample size of X and y. They should be the same!")
  } else if (ncol(y) != 1){
    stop("Error: Please check the dimension of y. It should be n*1 vector!")
  }
  # ---- check alpha0 range
  if (alpha0 <= 0 || alpha0 > 1){
    stop("Error: The range of alpha shoud be in (0,1]!")
  }
  # ---- check lambda0
  if (lambda0 <= 0){
    stop("Error: Lambda shoud be lager than 0!")
  }
  # ---- check gamma
  if (gamma <= 0){
    stop("Error: Gamma shoud be lager than 0!")
  }
  # --- check method
  if (method != "APGD" && method != "CVX"){
    stop("Error: The current methods must be 'APGD' or 'CVX'!")
  }

  # if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
  library(glmnet)
  ####### Setting:
  n <- nrow(X)
  p <- ncol(X)
  # M: shape parameter
  M <- IQR(t(X) %*% y, na.rm = TRUE)/1.345
  if (quiet == FALSE){
    print(paste("Start calculating beta using HuberENET by :", method))
  }
  # - APGD
  if (method == "APGD"){
    beta_est <- matrix(NA, nrow = ncol(X), ncol = niter+2)
    colnames(beta_est) <- paste("beta", seq(0, niter+1))
    obj <- numeric(niter+1)
    beta_est[,"beta 0"] <- beta_est[,"beta 1"] <- 0
    for (k in 1:niter){
      beta_k <- beta_est[,paste("beta", k)]
      z <- y - X %*% beta_k
      obj[k] <- sum(HuberFunction(z, M)) + lambda0*alpha0*sum(abs(beta_k)) + 0.5*lambda0*(1-alpha0)*t(beta_k)%*%beta_k
      ksi <- beta_est[,paste("beta", k)] + k/(k+3) * (beta_est[,paste("beta", k)] - beta_est[,paste("beta", k-1)])
      H_grad <- Huber_grid(y - X %*% ksi, M)
      g_grad_ksi <- - t(X) %*% H_grad + lambda0 * (1-alpha0) * ksi
      g_grad_ksi <- data.matrix(g_grad_ksi)
      while (TRUE){
        theta <- ksi - gamma * g_grad_ksi
        beta_prox <- sign(theta) * pmax(abs(theta)-gamma*lambda0*alpha0, 0)
        H_ksi <- HuberFunction(y - X %*% ksi, M)
        g_ksi <- sum(H_ksi) + 0.5*lambda0*(1-alpha0)*t(ksi)%*%ksi
        H_beta_prox <- HuberFunction(y - X %*% beta_prox, M)
        g_beta_prox <- sum(H_beta_prox) + 0.5*lambda0*(1-alpha0)*t(beta_prox)%*%beta_prox
        g_hat_upper <- g_ksi + t(g_grad_ksi)%*%(beta_prox-ksi) + 0.5/gamma*sum((beta_prox-ksi)^2)
        if (g_beta_prox > g_hat_upper){
          gamma <- 0.5*gamma
        } else {
          break
        }
      }
      beta_est[,paste("beta", k+1)] <- beta_prox
      z <- y - X %*% beta_prox
      obj[k+1] <-  sum(HuberFunction(z, M)) + lambda0*alpha0*sum(abs(beta_prox)) + 0.5*lambda0*(1-alpha0)*t(beta_prox)%*%beta_prox
      diff_beta <- beta_prox - beta_est[,paste("beta", k)]
      temp_beta <- (sum(abs(diff_beta) >= crit_beta) == 0)
      temp_obj <- (abs(obj[k+1] - obj[k]) < crit_obj)
      temp <- (1 - temp_beta) * (1 - temp_obj)
      if (k>1 && temp==0){
        break
      }
    }
    if (k == niter){
      print("Waining: Setting number of iterations doesn't reach to the convergency. Please set larger 'niter'!")
    }
    beta_hat <- beta_est[,k]
    if (quiet==FALSE && sum(beta_hat != 0) == 0){
      print("Warning: All estimated regression coefficients are 0s. Please check the size of lambda and input files!")
    }
  } else {   # - CVX
    # if (!requireNamespace("CVXR", quietly = TRUE)) install.packages("CVXR")
    library(CVXR)
    beta11 <- Variable(p)
    obj11 <- sum(huber(y - X%*%beta11, M)) + lambda0*alpha0*p_norm(beta11, 1) +
      0.5*lambda0*(1-alpha0)*sum(beta11^2)
    prob <- Problem(Minimize(obj11))
    result <- CVXR::solve(prob)
    beta_hat <- result[[1]]
  }

  return(beta_hat)
}
