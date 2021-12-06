\name{GraphicalModel}
\alias{GraphicalModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a covariance matrix from the specific graph (Adj) based on a Gaussian graphical model.
}
\description{
    %Simulate a covariance matrix from the specific graph (Adj) based on a Gaussian graphical model.
}
\usage{
GraphicalModel(Adj, a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Adj }{the adjacency matrix of network structure}
  \item{a1,a2,b1,b2 }{parameters for constructing domain [a1, a2] union [b1, b2]
#					default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{Sigma }{the covariance matrix of genes according to network structure.}

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

## Using "ConstructNetwork" to genarate network
N <- 200
Adj <- ConstructNetwork(N, "BAN")

sigma1 = GraphicalModel(Adj)
