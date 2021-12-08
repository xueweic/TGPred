\name{CalculateLaplacian}
\alias{CalculateLaplacian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Laplacian matrix and symmetric normalized Laplacian matrix from an adjacency matrix.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalculateLaplacian(Adj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Adj}{
    the adjacency matrix of network structure.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

    \item{L }{The Laplacian matrix for network structure.}
    \item{L_norm }{The symmetric normalized Laplacian matrix for network structure.}

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
  
  res <- CalculateLaplacian(Adj)
  
  L <- res$L
  
  L_norm <- res$L_norm
  
}
