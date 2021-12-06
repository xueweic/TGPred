\name{HuberLasso_SP}
\alias{HuberLasso_SP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate selection probability using HuberLasso function solving by APGD
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HuberLasso_SP(X, y, n_lambda, B = 500, gamma = 1000, niter = 2000, crit_beta = 1e-04, crit_obj = 1e-08, timer = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
expressional levels of n_genes target genes (TGs)
}
  \item{y}{
expressional levels of a transcription factor (TF)
}
  \item{n_lambda}{
the number of lambdas
}
  \item{B}{
the number of half-sample resamplings used to calculate selection probabilities of genes. default: 500
}
  \item{gamma}{
initial value of gamma in APGD. default: 1000
}
  \item{niter}{
the maximum number of APGD to solve Huber Lasso regression. default: 2000
}
  \item{crit_beta}{
converge criterion of change of beta. default: 1e-4
}
  \item{crit_obj}{
converge criterion of change of objective function. default: 1e-8
}
  \item{timer}{
decide if exist the output report. default: FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{SP }{n_genes length vector of selection probability.}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, y, n_lambda, B = 500, gamma = 1000, niter = 2000,
    crit_beta = 1e-04, crit_obj = 1e-08, timer = TRUE)
{
    X.ori <- data.matrix(X)
    X <- scale(X.ori)
    y.ori <- data.matrix(y)
    y <- scale(y.ori)
    if (nrow(X) != nrow(y)) {
        stop("Error: Please check the sample size of X and y. They should be the same!")
    }
    else if (ncol(y) != 1) {
        stop("Error: Please check the dimension of y. It should be n*1 vector!")
    }
    library(glmnet)
    n <- nrow(X)
    p <- ncol(X)
    M <- IQR(t(X) \%*\% y)/1.345
    print("Start calculating selection probability using HuberLasso by APGD method:")
    SP.LambdaAlpha <- matrix(NA, nrow = p, ncol = n_lambda)
    flag <- 0
    start.time <- proc.time()
    lambda_set <- Lambda_grid(X.ori, y.ori, n_lambda, alpha = 1,
        loss_func = "Huber")
    for (i.lambda in 1:length(lambda_set)) {
        lambda <- lambda_set[i.lambda]
        flag <- flag + 1
        print(paste0("Case ", flag, " : lambda ", lambda))
        beta_hat.LambdaAlpha <- matrix(NA, nrow = p, ncol = B)
        for (i.b in 1:B) {
            pos <- sample(1:n, n/2)
            y.sub <- y[pos, ]
            X.sub <- X[pos, ]
            beta_hat_APGD <- HuberLasso_Beta(X, y, lambda, method = "APGD",
                gamma = 1000, niter = 2000, crit_beta = 1e-04,
                crit_obj = 1e-08, quiet = TRUE)
            beta_hat.LambdaAlpha[, i.b] <- beta_hat_APGD
        }
        Proportion.LambdaAlpha <- (beta_hat.LambdaAlpha != 0) *
            1
        SP.LambdaAlpha[, flag] <- rowMeans(Proportion.LambdaAlpha)
    }
    SP.HuberLasso <- apply(SP.LambdaAlpha, 1, max)
    print("Done with HuberLasso by APGD method!")
    end.time <- proc.time()
    if (timer == TRUE) {
        diff <- end.time - start.time
        cost_h <- round(diff[3]/3600, 4)
        cost_min <- round(diff[3]/60, 4)
        cost_s <- round(diff[3], 4)
        print(paste("The calculate time:", cost_h, "h", cost_min,
            "min", cost_s, "s."))
    }
    return(SP.HuberLasso)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
