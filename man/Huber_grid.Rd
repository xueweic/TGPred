\name{Huber_grid}
\alias{Huber_grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Calculate the gridient of Huber function for an input value z.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Huber_grid(z, M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z }{Input value z.}
  \item{M }{Shape parameter, which is defaulted to be one-tenth of the interquartile range (IRQ).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 The gridient of Huber function for an input value z
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
function (delta, M)
{
    ifelse(abs(delta) < M, 2 * delta, 2 * M * sign(delta))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
