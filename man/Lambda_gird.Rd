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
  \item{alpha }{the proportion of l1 norm affects (the numerical values of nonzero coefficients). Its in range (0,1].}
  \item{loss_func }{either "Huber" or "MSE". If loss_func = "Huber", the loss function in penalized regression model is Huber function. If loss_func = "MSE", the loss function in penalized regression model is mean squared errors.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{lambdas }{n_lambda length vector of lambdas according to the alpha you provided.}

}
\examples{
alpha <- 0.5
    
n_lambda <- 10
    
lambda_set <- Lambda_grid(X, y, n_lambda, alpha, loss_func = "Huber")
    
lambda_set <- Lambda_grid(X, y, n_lambda, alpha, loss_func = "MSE")
}
