\name{ENET_Beta}
\alias{ENET_Beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Estimate beta_hat using Elastic Net function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ENET_Beta(X, y, lambda0, alpha0, method = "APGD", gamma = 1000, niter = 2000, crit_beta = 1e-04, crit_obj = 1e-08, quiet = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
expressional levels of n_genes target genes (TGs)
}
  \item{y}{
expressional levels of a transcription factor (TF)
}
  \item{lambda0}{
lambda0 one of parameters in Elastic Net regression, which controls the number of nonzero coefficients.
}
  \item{alpha0}{
alpha0 one of parameters in Elastic Net regression, which controls the numerical values of nonzero coefficients.
}
  \item{method}{
method The current methods must be 'APGD' or 'CVX'
}
  \item{gamma}{
gamma initial value of gamma in APGD. default: 1000
}
  \item{niter}{
niter the maximum number of APGD to solve Elastic Net regression. default: 2000
}
  \item{crit_beta}{
crit_beta converge criterion of change of beta. default: 1e-4
}
  \item{crit_obj}{
crit_obj converge criterion of change of objective function. default: 1e-8
}
  \item{quiet}{
quiet decide if exist the output report. default: FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{beta }{n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in ENET regression.}
}
\examples{

lambda0 = 200
alpha0 = 0.5
beta_hat_APGD <- ENET_Beta(X, y, lambda0, alpha0, method="APGD")
library("CVXR")
beta_hat_CVX <- ENET_Beta(X, y, lambda0, alpha0,method="CVX")
}


