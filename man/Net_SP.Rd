\name{Net_SP}
\alias{Net_SP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate selection probability using Net function solving by APGD
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Net_SP(X, y, Adj, alphas, n_lambda, B = 500, gamma = 1000, niter = 2000, crit_beta = 1e-04, crit_obj = 1e-08, timer = TRUE)
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
  \item{alphas}{
the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
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
the maximum number of APGD to solve Net regression. default: 2000
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

