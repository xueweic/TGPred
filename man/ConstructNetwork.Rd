\name{ConstructNetwork}
\alias{ConstructNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Construct the network structure from either Hierarchical Network or Barabasi-Albert Network in simulation studies.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ConstructNetwork(N_genes, Structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N_genes }{the number of genes}
  \item{structure }{"HN": Hierarchical Network or  "BAN": Barabasi-Albert Network}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
  \item{Adj }{the adjacency matrix of network structure}

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

N <- 200
Adj <- ConstructNetwork(N, "BAN")
Adj <- ConstructNetwork(N, "HN")
}
