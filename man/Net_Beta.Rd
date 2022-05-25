\name{Net_Beta}
\alias{Net_Beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate beta_hat using Net function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Net_Beta(X, y, Adj, lambda0, alpha0, method = "APGD", gamma = 1000, niter = 2000, crit_beta = 1e-04, crit_obj = 1e-08, quiet = FALSE, if.scale=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
expressional levels of n_genes target genes (TGs)
}
  \item{y}{
expressional levels of a transcription factor (TF)
}
  \item{Adj}{
the adjacency matrix of network structure.
}
  \item{lambda0}{
one of parameters in Net regression, which controls the number of nonzero coefficients.
}
  \item{alpha0}{
one of parameters in Net regression, which controls the numerical values of nonzero coefficients.
}
  \item{method}{
The current methods must be 'APGD' or 'CVX'
}
  \item{gamma}{
initial value of gamma in APGD. default: 1000
}
  \item{niter}{
the maximum number of APGD to solve HuberNet regression. default: 2000
}
  \item{crit_beta}{
converge criterion of change of beta. default: 1e-4
}
  \item{crit_obj}{
converge criterion of change of objective function. default: 1e-8
}
  \item{quiet}{
decide if exist the output report. default: FALSE
}
  \item{if.scale}{
if.scale decide if scale the expression levels. default: FALSE
} 
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{beta }{n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in Net regression.}

}
\examples{
lambda0 = 200
  
alpha0 = 0.5
  
beta_hat_APGD <- Net_Beta(X, y, Adj, lambda0, alpha0,method="APGD", if.scale=TRUE)
  
library("CVXR")
  
beta_hat_CVX <- Net_Beta(X, y, Adj,lambda0, alpha0,method="CVX", if.scale=TRUE)

}

