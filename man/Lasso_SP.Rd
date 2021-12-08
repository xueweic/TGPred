\name{Lasso_SP}
\alias{Lasso_SP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
stimate selection probability using Lasso function solving by APGD
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lasso_SP(X, y, n_lambda, B = 500, gamma = 1000, niter = 2000, crit_beta = 1e-04, crit_obj = 1e-08, timer = TRUE)
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
the maximum number of APGD to solve Lasso regression. default: 2000
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
\examples{
n_lambda <- 50

SP_Lasso = Lasso_SP(X, y, n_lambda, B=B0, gamma=1000, niter=2000, timer=FALSE)
}

