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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
