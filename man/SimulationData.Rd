\name{SimulationData}
\alias{SimulationData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Simulate y and X from a given network structure.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SimulationData(N_samples, N_genes, Adj, Sigma, method, beta0=NULL, beta_true=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N_samples }{the number of sample size.}
  \item{N_genes }{the number of traget genes.}
  \item{Adj }{the adjacency matrix of network structure. Adjacency matrix must be a N_genes*N_genes dimensional symmetric matrix, the elements equal 1 indicates two genes are connected. If you consider Barabasi-Albert Network or Hierarchical Network in the article, you can directly use "ConstructNetwork" function to get the adjacency matrix.}
  \item{Sigma }{the covariance matrix of target genes according to network structure. You can directly use "GraphicalModel" function to get the covariance matrix.}
  \item{method }{"HN": by Hierarchical Network, "BAN": by Barabasi-Albert Network or "DIY": by user designed}
  \item{beta0 }{numeric value of effect size in simulation settings.
#		   default: NULL; if method is "HN" or "BAN", input a nunerical value.}
  \item{beta_true }{numeric matrix with the dimension of N_genes * 1 in simulation settings.
#		   default: NULL; if method is "DIY", input a nunerical matrix (N_genes * 1).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

 \item{y }{expression levels of a transcription factor (TF)}
 \item{X }{expression levels of n_genes target genes (TGs)}
 \item{beta }{true regulated effect beta for N_genes TGs.}

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

N.sample <- 300
N <- 200
Adj = ConstructNetwork(N, "BAN")
sigma1 = GraphicalModel(Adj)


res = SimulationData(N.sample,N,Adj,sigma1,"BAN", beta0 = 10)
y = res$y
X = res$X
beta1 = res$beta


