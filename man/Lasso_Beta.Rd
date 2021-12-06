\name{Lasso_Beta}
\alias{Lasso_Beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate beta_hat using Lasso function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lasso_Beta(X, y, lambda0, method = "APGD", gamma = 1000, niter = 2000, crit_beta = 1e-04, crit_obj = 1e-08, quiet = FALSE)
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
one of parameters in Lasso regression, which controls the number of nonzero coefficients.
}
  \item{method}{
The current methods must be 'APGD' or 'CVX'
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
  \item{quiet}{
decide if exist the output report. default: FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{beta }{n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in Lasso regression.}

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
