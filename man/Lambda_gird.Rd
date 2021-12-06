\name{Lambda_gird}
\alias{Lambda_gird}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Simulate a grid set of lambdas for a given alpha in penalized regression.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lambda_gird(X, y, n_lambda, alpha, loss_func)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{expression levels of n_genes target genes (TGs).}
  \item{y }{expression levels of a transcription factor (TF).}
  \item{n_lambda }{the number of lambdas.}
  \item{alpha }{the proportion of l1 norm affects (the numerical values of nonzero coefficients). It's in range (0,1].}
  \item{loss_func }{either "Huber" or "MSE".
  # If loss_func = "Huber", the loss function in penalized regression model is Huber function.
  # If loss_func = "MSE", the loss function in penalized regression model is mean squared errors.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{lambdas }{n_lambda length vector of lambdas according to the alpha you provided.}

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
